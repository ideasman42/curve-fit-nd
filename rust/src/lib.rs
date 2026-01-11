//! Curve Fitting Library
//!
//! Attempt to fit a cubic bezier curve to a set of points.
//! This is a Rust port of the C curve-fit-nd library.

use libc::{c_double, c_float, c_uint};
use std::ptr;

mod math_vector;
mod min_heap;

mod curve_fit_cubic;
mod curve_fit_cubic_refit;
mod curve_fit_from_polys;

pub use curve_fit_from_polys::{fit_poly_single, fit_poly_single_with_corners, TraceMode};

/// Calculation flags matching the C API.
pub const CURVE_FIT_CALC_HIGH_QUALITY: c_uint = 1 << 0;
pub const CURVE_FIT_CALC_CYCLIC: c_uint = 1 << 1;

/// Fit cubic bezier curves to points using the refit algorithm.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_refit_db(
    points: *const c_double,
    points_len: c_uint,
    dims: c_uint,
    error_threshold: c_double,
    calc_flag: c_uint,
    corners: *const c_uint,
    corners_len: c_uint,
    corner_angle: c_double,
    r_cubic_array: *mut *mut c_double,
    r_cubic_array_len: *mut c_uint,
    r_cubic_orig_index: *mut *mut c_uint,
    r_corner_index_array: *mut *mut c_uint,
    r_corner_index_len: *mut c_uint,
) -> i32 {
    // Validate inputs
    if points.is_null() || points_len < 2 || dims < 2 {
        return -1;
    }

    let dims_usize = dims as usize;
    let is_cyclic = (calc_flag & CURVE_FIT_CALC_CYCLIC) != 0;
    let use_high_quality = (calc_flag & CURVE_FIT_CALC_HIGH_QUALITY) != 0;

    // Create a slice view directly from the C pointer (zero-copy).
    let points_slice: &[f64] = std::slice::from_raw_parts(
        points,
        (points_len as usize) * dims_usize,
    );

    // Convert corners to a Rust slice (if provided)
    let corners_slice: Option<Vec<usize>> = if !corners.is_null() && corners_len > 0 {
        let raw_corners = std::slice::from_raw_parts(corners, corners_len as usize);
        Some(raw_corners.iter().map(|&c| c as usize).collect())
    } else {
        None
    };

    // Call the Rust implementation
    let (result, orig_indices) = fit_poly_single_with_corners(
        points_slice,
        dims_usize,
        is_cyclic,
        error_threshold,
        corner_angle,
        use_high_quality,
        corners_slice.as_deref(),
    );

    // result is a flat array: [h_in_0, p_0, h_out_0, h_in_1, p_1, h_out_1, ...]
    // Each of h_in, p, h_out has `dims` components
    // So result.len() = num_knots * 3 * dims
    let num_knots = orig_indices.len();
    if num_knots == 0 {
        *r_cubic_array = ptr::null_mut();
        *r_cubic_array_len = 0;
        if !r_cubic_orig_index.is_null() {
            *r_cubic_orig_index = ptr::null_mut();
        }
        if !r_corner_index_array.is_null() {
            *r_corner_index_array = ptr::null_mut();
        }
        if !r_corner_index_len.is_null() {
            *r_corner_index_len = 0;
        }
        return 0;
    }

    // Allocate output arrays using libc malloc (so C code can free them)
    let cubic_data_len = result.len(); // num_knots * 3 * dims
    let cubic_array_ptr = libc::malloc(cubic_data_len * std::mem::size_of::<c_double>()) as *mut c_double;
    if cubic_array_ptr.is_null() {
        return -1;
    }

    // Copy cubic data
    let cubic_slice = std::slice::from_raw_parts_mut(cubic_array_ptr, cubic_data_len);
    cubic_slice.copy_from_slice(&result);

    *r_cubic_array = cubic_array_ptr;
    *r_cubic_array_len = num_knots as c_uint;

    // Allocate and fill original index array
    if !r_cubic_orig_index.is_null() {
        let orig_index_ptr = libc::malloc(num_knots * std::mem::size_of::<c_uint>()) as *mut c_uint;
        if orig_index_ptr.is_null() {
            libc::free(cubic_array_ptr as *mut libc::c_void);
            return -1;
        }
        let orig_index_slice = std::slice::from_raw_parts_mut(orig_index_ptr, num_knots);
        for (i, idx) in orig_index_slice.iter_mut().enumerate() {
            *idx = orig_indices[i] as c_uint;
        }
        *r_cubic_orig_index = orig_index_ptr;
    }

    // Corner indices (not yet fully implemented)
    if !r_corner_index_array.is_null() {
        *r_corner_index_array = ptr::null_mut();
    }
    if !r_corner_index_len.is_null() {
        *r_corner_index_len = 0;
    }

    0
}

/// Float version of #curve_fit_cubic_to_points_refit.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_refit_fl(
    points: *const c_float,
    points_len: c_uint,
    dims: c_uint,
    error_threshold: c_float,
    calc_flag: c_uint,
    corners: *const c_uint,
    corners_len: c_uint,
    corner_angle: c_float,
    r_cubic_array: *mut *mut c_float,
    r_cubic_array_len: *mut c_uint,
    r_cubic_orig_index: *mut *mut c_uint,
    r_corner_index_array: *mut *mut c_uint,
    r_corner_index_len: *mut c_uint,
) -> i32 {
    // Convert float inputs to double, call double version, convert back
    if points.is_null() || points_len < 2 || dims < 2 {
        return -1;
    }

    let dims_usize = dims as usize;
    let points_slice = std::slice::from_raw_parts(points, (points_len as usize) * dims_usize);
    let points_double: Vec<c_double> = points_slice.iter().map(|&x| x as c_double).collect();

    let mut cubic_array_db: *mut c_double = ptr::null_mut();
    let mut cubic_array_len: c_uint = 0;
    let mut cubic_orig_index: *mut c_uint = ptr::null_mut();
    let mut corner_index_array: *mut c_uint = ptr::null_mut();
    let mut corner_index_len: c_uint = 0;

    let result = curve_fit_cubic_to_points_refit_db(
        points_double.as_ptr(),
        points_len,
        dims,
        error_threshold as c_double,
        calc_flag,
        corners,
        corners_len,
        corner_angle as c_double,
        &mut cubic_array_db,
        &mut cubic_array_len,
        if r_cubic_orig_index.is_null() { ptr::null_mut() } else { &mut cubic_orig_index },
        if r_corner_index_array.is_null() { ptr::null_mut() } else { &mut corner_index_array },
        if r_corner_index_len.is_null() { ptr::null_mut() } else { &mut corner_index_len },
    );

    if result != 0 {
        return result;
    }

    // Convert double results back to float
    if !cubic_array_db.is_null() && cubic_array_len > 0 {
        let cubic_data_len = (cubic_array_len as usize) * 3 * dims_usize;
        let cubic_array_fl = libc::malloc(cubic_data_len * std::mem::size_of::<c_float>()) as *mut c_float;
        if cubic_array_fl.is_null() {
            libc::free(cubic_array_db as *mut libc::c_void);
            if !cubic_orig_index.is_null() {
                libc::free(cubic_orig_index as *mut libc::c_void);
            }
            return -1;
        }

        let db_slice = std::slice::from_raw_parts(cubic_array_db, cubic_data_len);
        let fl_slice = std::slice::from_raw_parts_mut(cubic_array_fl, cubic_data_len);
        for (i, &val) in db_slice.iter().enumerate() {
            fl_slice[i] = val as c_float;
        }

        libc::free(cubic_array_db as *mut libc::c_void);
        *r_cubic_array = cubic_array_fl;
    } else {
        *r_cubic_array = ptr::null_mut();
    }

    *r_cubic_array_len = cubic_array_len;

    if !r_cubic_orig_index.is_null() {
        *r_cubic_orig_index = cubic_orig_index;
    }
    if !r_corner_index_array.is_null() {
        *r_corner_index_array = corner_index_array;
    }
    if !r_corner_index_len.is_null() {
        *r_corner_index_len = corner_index_len;
    }

    0
}

// ============================================================================
// Non-refit curve fitting (delegates to refit version)
// ============================================================================

/// Fit cubic bezier curves to points (non-refit version).
///
/// This delegates to the refit version with corner_angle set to PI (no corner detection).
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_db(
    points: *const c_double,
    points_len: c_uint,
    dims: c_uint,
    error_threshold: c_double,
    calc_flag: c_uint,
    corners: *const c_uint,
    corners_len: c_uint,
    r_cubic_array: *mut *mut c_double,
    r_cubic_array_len: *mut c_uint,
    r_cubic_orig_index: *mut *mut c_uint,
    r_corner_index_array: *mut *mut c_uint,
    r_corner_index_len: *mut c_uint,
) -> i32 {
    // Delegate to refit version with corner_angle = PI (no corner detection)
    curve_fit_cubic_to_points_refit_db(
        points,
        points_len,
        dims,
        error_threshold,
        calc_flag,
        corners,
        corners_len,
        std::f64::consts::PI,
        r_cubic_array,
        r_cubic_array_len,
        r_cubic_orig_index,
        r_corner_index_array,
        r_corner_index_len,
    )
}

/// Float version of #curve_fit_cubic_to_points.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_fl(
    points: *const c_float,
    points_len: c_uint,
    dims: c_uint,
    error_threshold: c_float,
    calc_flag: c_uint,
    corners: *const c_uint,
    corners_len: c_uint,
    r_cubic_array: *mut *mut c_float,
    r_cubic_array_len: *mut c_uint,
    r_cubic_orig_index: *mut *mut c_uint,
    r_corner_index_array: *mut *mut c_uint,
    r_corner_index_len: *mut c_uint,
) -> i32 {
    // Delegate to refit version with corner_angle = PI (no corner detection)
    curve_fit_cubic_to_points_refit_fl(
        points,
        points_len,
        dims,
        error_threshold,
        calc_flag,
        corners,
        corners_len,
        std::f32::consts::PI,
        r_cubic_array,
        r_cubic_array_len,
        r_cubic_orig_index,
        r_corner_index_array,
        r_corner_index_len,
    )
}

// ============================================================================
// Single segment curve fitting
// ============================================================================

/// Fit a single cubic bezier curve to a set of points.
///
/// Returns 0 if threshold was met (success), 1 if not met, -1 on error.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_single_db(
    points: *const c_double,
    points_len: c_uint,
    points_length_cache: *const c_double,
    dims: c_uint,
    error_threshold: c_double,
    tan_l: *const c_double,
    tan_r: *const c_double,
    r_handle_l: *mut c_double,
    r_handle_r: *mut c_double,
    r_error_sq: *mut c_double,
    r_error_index: *mut c_uint,
) -> i32 {
    // Validate inputs
    if points.is_null() || tan_l.is_null() || tan_r.is_null() || dims < 2 {
        return -1;
    }
    if points_len < 2 {
        return -1;
    }

    let dims_usize = dims as usize;
    let points_len_usize = points_len as usize;
    let error_threshold_sq = error_threshold * error_threshold;

    // Create slice views
    let points_slice: &[f64] = std::slice::from_raw_parts(
        points,
        points_len_usize * dims_usize,
    );
    let tan_l_slice: &[f64] = std::slice::from_raw_parts(tan_l, dims_usize);
    let tan_r_slice: &[f64] = std::slice::from_raw_parts(tan_r, dims_usize);

    // Handle optional length cache
    let length_cache_vec: Vec<f64>;
    let length_cache: &[f64] = if points_length_cache.is_null() {
        // Calculate length cache if not provided
        length_cache_vec = calculate_length_cache(points_slice, dims_usize);
        &length_cache_vec
    } else {
        std::slice::from_raw_parts(points_length_cache, points_len_usize)
    };

    // Call the Rust implementation
    let ((error_sq, error_index), handle_l, handle_r, threshold_met) =
        curve_fit_cubic::curve_fit_cubic_to_points_single(
            points_slice,
            dims_usize,
            length_cache,
            tan_l_slice,
            tan_r_slice,
            error_threshold_sq,
        );

    // Write results
    if !r_handle_l.is_null() {
        let handle_l_out = std::slice::from_raw_parts_mut(r_handle_l, dims_usize);
        handle_l_out.copy_from_slice(&handle_l[..dims_usize]);
    }
    if !r_handle_r.is_null() {
        let handle_r_out = std::slice::from_raw_parts_mut(r_handle_r, dims_usize);
        handle_r_out.copy_from_slice(&handle_r[..dims_usize]);
    }
    if !r_error_sq.is_null() {
        *r_error_sq = error_sq;
    }
    if !r_error_index.is_null() {
        *r_error_index = error_index as c_uint;
    }

    // Return 0 for threshold met (success), 1 for not met.
    if threshold_met { 0 } else { 1 }
}

/// Float version of #curve_fit_cubic_to_points_single.
///
/// Returns 0 if threshold was met (success), 1 if not met, -1 on error.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_cubic_to_points_single_fl(
    points: *const c_float,
    points_len: c_uint,
    points_length_cache: *const c_float,
    dims: c_uint,
    error_threshold: c_float,
    tan_l: *const c_float,
    tan_r: *const c_float,
    r_handle_l: *mut c_float,
    r_handle_r: *mut c_float,
    r_error_sq: *mut c_float,
    r_error_index: *mut c_uint,
) -> i32 {
    if points.is_null() || tan_l.is_null() || tan_r.is_null() || dims < 2 {
        return -1;
    }

    let dims_usize = dims as usize;
    let points_len_usize = points_len as usize;

    // Convert to double
    let points_slice = std::slice::from_raw_parts(points, points_len_usize * dims_usize);
    let points_double: Vec<c_double> = points_slice.iter().map(|&x| x as c_double).collect();

    let length_cache_double: Vec<c_double> = if points_length_cache.is_null() {
        Vec::new()
    } else {
        let cache_slice = std::slice::from_raw_parts(points_length_cache, points_len_usize);
        cache_slice.iter().map(|&x| x as c_double).collect()
    };

    let tan_l_slice = std::slice::from_raw_parts(tan_l, dims_usize);
    let tan_l_double: Vec<c_double> = tan_l_slice.iter().map(|&x| x as c_double).collect();
    let tan_r_slice = std::slice::from_raw_parts(tan_r, dims_usize);
    let tan_r_double: Vec<c_double> = tan_r_slice.iter().map(|&x| x as c_double).collect();

    let mut handle_l_double: Vec<c_double> = vec![0.0; dims_usize];
    let mut handle_r_double: Vec<c_double> = vec![0.0; dims_usize];
    let mut error_sq_double: c_double = 0.0;
    let mut error_index: c_uint = 0;

    let result = curve_fit_cubic_to_points_single_db(
        points_double.as_ptr(),
        points_len,
        if length_cache_double.is_empty() { ptr::null() } else { length_cache_double.as_ptr() },
        dims,
        error_threshold as c_double,
        tan_l_double.as_ptr(),
        tan_r_double.as_ptr(),
        handle_l_double.as_mut_ptr(),
        handle_r_double.as_mut_ptr(),
        &mut error_sq_double,
        &mut error_index,
    );

    // Return early only on error, not on threshold not met.
    if result < 0 {
        return result;
    }

    // Convert back to float
    if !r_handle_l.is_null() {
        let handle_l_out = std::slice::from_raw_parts_mut(r_handle_l, dims_usize);
        for (i, &val) in handle_l_double.iter().enumerate() {
            handle_l_out[i] = val as c_float;
        }
    }
    if !r_handle_r.is_null() {
        let handle_r_out = std::slice::from_raw_parts_mut(r_handle_r, dims_usize);
        for (i, &val) in handle_r_double.iter().enumerate() {
            handle_r_out[i] = val as c_float;
        }
    }
    if !r_error_sq.is_null() {
        *r_error_sq = error_sq_double as c_float;
    }
    if !r_error_index.is_null() {
        *r_error_index = error_index;
    }

    result
}

// ============================================================================
// Corner detection
// ============================================================================

/// Detect corners in a point sequence.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_corners_detect_db(
    _points: *const c_double,
    _points_len: c_uint,
    _dims: c_uint,
    _radius_min: c_double,
    _radius_max: c_double,
    _samples_max: c_uint,
    _angle_threshold: c_double,
    r_corners: *mut *mut c_uint,
    r_corners_len: *mut c_uint,
) -> i32 {
    // TODO: Implement corner detection
    // For now, return empty corner array
    if !r_corners.is_null() {
        *r_corners = ptr::null_mut();
    }
    if !r_corners_len.is_null() {
        *r_corners_len = 0;
    }
    0
}

/// Float version of #curve_fit_corners_detect.
///
/// # Safety
/// This function is meant to be called from C code. All pointer arguments must be valid.
#[no_mangle]
pub unsafe extern "C" fn curve_fit_corners_detect_fl(
    _points: *const c_float,
    _points_len: c_uint,
    _dims: c_uint,
    _radius_min: c_float,
    _radius_max: c_float,
    _samples_max: c_uint,
    _angle_threshold: c_float,
    r_corners: *mut *mut c_uint,
    r_corners_len: *mut c_uint,
) -> i32 {
    // TODO: Implement corner detection
    // For now, return empty corner array
    if !r_corners.is_null() {
        *r_corners = ptr::null_mut();
    }
    if !r_corners_len.is_null() {
        *r_corners_len = 0;
    }
    0
}

// ============================================================================
// Helper functions
// ============================================================================

/// Calculate the length cache (distances between consecutive points).
fn calculate_length_cache(points: &[f64], dims: usize) -> Vec<f64> {
    let points_len = points.len() / dims;
    let mut cache = Vec::with_capacity(points_len);
    cache.push(0.0); // First point has no previous distance

    for i in 1..points_len {
        let mut dist_sq = 0.0;
        for d in 0..dims {
            let diff = points[i * dims + d] - points[(i - 1) * dims + d];
            dist_sq += diff * diff;
        }
        cache.push(dist_sq.sqrt());
    }

    cache
}
