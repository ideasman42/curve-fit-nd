
///
/// Perform cubic curve fitting
///
/// This module takes a complete polygon and optimizes curve fitting
/// and optionally corner calculation,
/// outputting a bezier curve that fits within an error margin.
///

/// Enable the refit pass for improved curve quality.
const USE_REFIT: bool = true;
/// Allow removing knots during the refit pass.
const USE_REFIT_REMOVE: bool = true;
/// Scale factor for corner detection error threshold.
/// The C implementation uses 3.0 for the collapse threshold, allowing corners
/// to be collapsed with up to 9x the squared error tolerance.
const CORNER_SCALE: f64 = 3.0;

use crate::math_vector::{
    add_vnvn,
    copy_vnvn,
    madd_vnvn_fl,
    normalize_vn,
    normalized_vnvn_with_len,
    sq,
    zero_vn,
    is_finite_vn,
};

use crate::min_heap;

// Import refit module types and functions
use super::curve_fit_cubic_refit::{
    INVALID,
    Knot,
    PointData,
    refine_remove,
    refine_refit,
    refine_corner,
};

/// Tracing mode for polygon extraction from raster images.
#[derive(Copy, Clone, PartialEq)]
pub enum TraceMode {
    /// Trace the outline (boundary) of shapes.
    Outline,
    /// Trace the centerline (skeleton) of shapes.
    Centerline,
}

/// Get a point slice from a flat points array.
#[inline]
fn get_point(points: &[f64], dims: usize, index: usize) -> &[f64] {
    let start = index * dims;
    &points[start..start + dims]
}

/// Get a mutable tangent slice from a flat tangents array.
#[inline]
fn get_tangent_mut(tangents: &mut [f64], dims: usize, tan_index: usize) -> &mut [f64] {
    let start = tan_index * dims;
    &mut tangents[start..start + dims]
}


/// Fit cubic bezier curves to a single polygon.
///
/// Takes a list of points (as a flat array) and fits bezier curves that approximate
/// the polygon within the specified error threshold.
///
/// # Arguments
/// * `points_orig` - The input points as a flat array [x0, y0, z0, x1, y1, z1, ...]
/// * `dims` - Number of dimensions per point (e.g., 2 for 2D, 3 for 3D)
/// * `is_cyclic` - Whether the curve is closed/cyclic
/// * `error_threshold` - Maximum allowed error for curve fitting
/// * `corner_angle` - Angle threshold for automatic corner detection (PI = no detection)
/// * `use_optimize_exhaustive` - Use exhaustive optimization (slower but better)
///
/// Returns a tuple of:
/// - A flat array of cubic bezier segments [h_in_x, h_in_y, ..., p_x, p_y, ..., h_out_x, h_out_y, ...]
/// - A list of original point indices corresponding to each bezier segment.
pub fn fit_poly_single(
    points_orig: &[f64],
    dims: usize,
    is_cyclic: bool,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
) -> (Vec<f64>, Vec<usize>) {
    fit_poly_single_with_corners(
        points_orig,
        dims,
        is_cyclic,
        error_threshold,
        corner_angle,
        use_optimize_exhaustive,
        None,
    )
}

/// Fit cubic bezier curves to a single polygon with optional pre-defined corners.
///
/// See [`fit_poly_single`] for basic documentation.
pub fn fit_poly_single_with_corners(
    points_orig: &[f64],
    dims: usize,
    is_cyclic: bool,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
    corners: Option<&[usize]>,
) -> (Vec<f64>, Vec<usize>) {
    use std::borrow::Cow;

    let points_len = points_orig.len() / dims;
    let knots_len = points_len;

    // For cyclic curves, duplicate the points array to allow contiguous slicing
    // across the start/end boundary (matching C behavior).
    // For non-cyclic curves, use the original slice directly (zero-copy).
    let points: Cow<'_, [f64]> = if is_cyclic {
        Cow::Owned([points_orig, points_orig].concat())
    } else {
        Cow::Borrowed(points_orig)
    };

    let mut knots: Vec<Knot> = Vec::with_capacity(knots_len);
    let mut knots_handle: Vec<min_heap::NodeHandle> =
        vec![min_heap::NodeHandle::INVALID; knots_len];

    let use_corner = corner_angle < ::std::f64::consts::PI;

    for i in 0..knots_len {
        assert!(is_finite_vn(get_point(points_orig, dims, i)));
        knots.push(Knot {
            next: i.wrapping_add(1),
            prev: i.wrapping_sub(1),
            index: i,
            no_remove: false,
            is_remove: false,
            is_corner: false,
            handles: [-1.0, -1.0], // dummy
            fit_error_sq_next: 0.0,
            tan: [i * 2, i * 2 + 1],
        });
    }

    if is_cyclic {
        let i_last = knots.len() - 1;
        knots[0].prev = i_last;
        knots[i_last].next = 0;
    } else {
        let i_last = knots.len() - 1;
        knots[0].prev = INVALID;
        knots[i_last].next = INVALID;

        knots[0].no_remove = true;
        knots[i_last].no_remove = true;
    }

    // All values will be written to, simplest to initialize to dummy values for now.
    // Double the cache for cyclic curves to match C behavior.
    let mut points_length_cache: Vec<f64> = vec![-1.0; points_len * if is_cyclic { 2 } else { 1 }];
    // Tangents: 2 per knot (incoming and outgoing), each with `dims` components
    let mut tangents: Vec<f64> = vec![-1.0; knots_len * 2 * dims];

    // Initialize tangents,
    // also set the values for knot handles since some may not collapse.

    if knots_len < 2 {
        for (i, k) in (&mut knots).iter_mut().enumerate() {
            zero_vn(get_tangent_mut(&mut tangents, dims, k.tan[0]));
            zero_vn(get_tangent_mut(&mut tangents, dims, k.tan[1]));
            k.handles[0] = 0.0;
            k.handles[1] = 0.0;
            points_length_cache[i] = 0.0;
        }
    } else if is_cyclic {
        let (mut tan_prev, mut len_prev) = normalized_vnvn_with_len(
            get_point(&points, dims, knots_len - 2),
            get_point(&points, dims, knots_len - 1),
            dims);

        let mut i_curr = knots.len() - 1;
        for i_next in 0..knots.len() {
            let k = &mut knots[i_curr];

            let (tan_next, len_next) = normalized_vnvn_with_len(
                get_point(&points, dims, i_curr),
                get_point(&points, dims, i_next),
                dims);
            points_length_cache[i_next] = len_next;

            let mut t = add_vnvn(&tan_prev[..dims], &tan_next[..dims], dims);
            normalize_vn(&mut t[..dims]);
            assert!(is_finite_vn(&t[..dims]));
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[0]), &t[..dims]);
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[1]), &t[..dims]);

            k.handles[0] = len_prev /  3.0;
            k.handles[1] = len_next / -3.0;

            tan_prev = tan_next;
            len_prev = len_next;
            i_curr = i_next;
        }
    } else {
        points_length_cache[0] = 0.0;
        let (mut tan_prev, mut len_prev) = normalized_vnvn_with_len(
            get_point(&points, dims, 0),
            get_point(&points, dims, 1),
            dims);
        points_length_cache[1] = len_prev;

        copy_vnvn(get_tangent_mut(&mut tangents, dims, knots[0].tan[0]), &tan_prev[..dims]);
        copy_vnvn(get_tangent_mut(&mut tangents, dims, knots[0].tan[1]), &tan_prev[..dims]);
        knots[0].handles[0] = len_prev /  3.0;
        knots[0].handles[1] = len_prev / -3.0;

        let mut i_curr = 1;
        for i_next in 2..knots.len() {
            let k = &mut knots[i_curr];
            let (tan_next, len_next) = normalized_vnvn_with_len(
                get_point(&points, dims, i_curr),
                get_point(&points, dims, i_next),
                dims);
            points_length_cache[i_next] = len_next;

            let mut t = add_vnvn(&tan_prev[..dims], &tan_next[..dims], dims);
            normalize_vn(&mut t[..dims]);
            assert!(is_finite_vn(&t[..dims]));
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[0]), &t[..dims]);
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[1]), &t[..dims]);

            k.handles[0] = len_prev /  3.0;
            k.handles[1] = len_next / -3.0;

            tan_prev = tan_next;
            len_prev = len_next;
            i_curr = i_next;
        }
        // use prev as next since they're copied above
        copy_vnvn(get_tangent_mut(&mut tangents, dims, knots[knots_len - 1].tan[0]), &tan_prev[..dims]);
        copy_vnvn(get_tangent_mut(&mut tangents, dims, knots[knots_len - 1].tan[1]), &tan_prev[..dims]);

        knots[knots_len - 1].handles[0] = len_prev /  3.0;
        knots[knots_len - 1].handles[1] = len_prev / -3.0;
    }

    // Duplicate the length cache for cyclic curves (matching C behavior).
    if is_cyclic {
        for i in 0..points_len {
            points_length_cache[i + points_len] = points_length_cache[i];
        }
    }

    // Initialize pre-defined corners and their tangents.
    // This overwrites the smooth tangents for corner points with separate
    // tangents pointing to adjacent points.
    if let Some(corner_indices) = corners {
        let knots_end = knots_len - 1;
        // For non-cyclic curves, skip first and last corners (handled as endpoints).
        let corners_start = if is_cyclic { 0 } else { 1 };
        let corners_len_clamped = if is_cyclic { corner_indices.len() } else { corner_indices.len().saturating_sub(1) };

        for corner_i in corners_start..corners_len_clamped {
            let i_curr = corner_indices[corner_i];
            let i_prev = if is_cyclic && i_curr == 0 { knots_end } else { i_curr - 1 };
            let i_next = if is_cyclic && i_curr == knots_end { 0 } else { i_curr + 1 };

            let k = &mut knots[i_curr];
            // Tangent towards previous point (for incoming handle).
            let (tan_prev, len_prev) = normalized_vnvn_with_len(
                get_point(&points, dims, i_prev),
                get_point(&points, dims, i_curr),
                dims);
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[0]), &tan_prev[..dims]);
            k.handles[0] = len_prev / 3.0;

            // Tangent towards next point (for outgoing handle).
            let (tan_next, len_next) = normalized_vnvn_with_len(
                get_point(&points, dims, i_curr),
                get_point(&points, dims, i_next),
                dims);
            copy_vnvn(get_tangent_mut(&mut tangents, dims, k.tan[1]), &tan_next[..dims]);
            k.handles[1] = len_next / -3.0;

            k.is_corner = true;
        }
    }

    let mut knots_len_remaining = knots.len();
    let pd = PointData {
        points: &*points,
        dims: dims,
        points_len: points_len,
        points_length_cache: &points_length_cache,
        tangents: &tangents,
    };

    // `curve_incremental_simplify_refit` can be called here, but it's very slow,
    // just remove all within the threshold first.
    refine_remove::curve_incremental_simplify(
        &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
        sq(error_threshold));

    if use_corner {
        refine_corner::curve_incremental_simplify_corners(
            &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
            sq(error_threshold), sq(error_threshold * CORNER_SCALE),
            corner_angle,
            );
    }

    debug_assert!(knots_len_remaining >= 2);

    if USE_REFIT {
        refine_refit::curve_incremental_simplify_refit(
            &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
            sq(error_threshold), use_optimize_exhaustive, USE_REFIT_REMOVE);
    }

    debug_assert!(knots_len_remaining >= 2);

    // Correct unused handle endpoints - not essential, but nice behavior.
    // For non-cyclic curves, the first knot's incoming handle and the last knot's
    // outgoing handle are unused. Set them to mirror their used counterparts.
    if !is_cyclic {
        // Find first active knot
        let k_first_index: usize = knots.iter().position(|k| !k.is_remove).unwrap();
        // Find last active knot
        let mut k_last_index = k_first_index;
        let mut k_idx = k_first_index;
        for _ in 0..knots_len_remaining {
            k_last_index = k_idx;
            k_idx = knots[k_idx].next;
        }
        knots[k_first_index].handles[0] = -knots[k_first_index].handles[1];
        knots[k_last_index].handles[1] = -knots[k_last_index].handles[0];
    }

    // Output: flat array of cubic bezier segments
    // Each segment has 3 points (handle_in, point, handle_out), each with `dims` components
    // Total size: knots_len_remaining * 3 * dims
    let mut cubic_array: Vec<f64> = Vec::with_capacity(knots_len_remaining * 3 * dims);
    let mut orig_indices: Vec<usize> = Vec::with_capacity(knots_len_remaining);

    /// Get a tangent slice from a flat tangents array.
    #[inline]
    fn get_tangent(tangents: &[f64], dims: usize, tan_index: usize) -> &[f64] {
        let start = tan_index * dims;
        &tangents[start..start + dims]
    }

    {
        let k_first_index: usize = {
            let mut i_search = INVALID;
            for (i, k) in knots.iter().enumerate() {
                if k.is_remove == false {
                    i_search = i;
                    break;
                }
            }
            debug_assert!(i_search != INVALID);
            i_search
        };

        let mut k_index = k_first_index;
        for _ in 0..knots_len_remaining {
            let k = &knots[k_index];
            let p = get_point(&points, dims, k.index);

            // handle_in = p + tangent[0] * handles[0]
            let handle_in = madd_vnvn_fl(p, get_tangent(&tangents, dims, k.tan[0]), k.handles[0], dims);
            cubic_array.extend_from_slice(&handle_in[..dims]);

            // point
            cubic_array.extend_from_slice(p);

            // handle_out = p + tangent[1] * handles[1]
            let handle_out = madd_vnvn_fl(p, get_tangent(&tangents, dims, k.tan[1]), k.handles[1], dims);
            cubic_array.extend_from_slice(&handle_out[..dims]);

            orig_indices.push(k.index);

            k_index = k.next;
        }
    }

    (cubic_array, orig_indices)
}
