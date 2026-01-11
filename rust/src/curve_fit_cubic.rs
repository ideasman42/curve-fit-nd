///
/// Curve Fitting - Single Segment
/// ===============================
///
/// This module fits a single cubic bezier curve to a set of points.
///
/// It implements multiple fitting strategies and selects the best result:
/// - Fallback: Simple 1/3 endpoint distance handles (baseline)
/// - Circular: Arc-based approximation for curved segments
/// - Offset: Perpendicular distance method for symmetric curves
/// - Least-squares: Optimal Bernstein polynomial fitting
///
/// The least-squares solution uses Newton-Raphson iteration to refine
/// the parameterization, and applies handle clamping to prevent
/// extreme control point positions.
///

use crate::math_vector::{
    len_squared_vnvn,
    len_vnvn,
    sub_vnvn,
    dot_vnvn,
    sq,
    DIMS_MAX,
};

// ============================================================================
// Cubic Type
// ============================================================================

/// A cubic bezier curve defined by four control points.
#[derive(Clone, Copy)]
struct Cubic {
    /// Start point (on the curve).
    p0: [f64; DIMS_MAX],
    /// First control point (handle from p0).
    p1: [f64; DIMS_MAX],
    /// Second control point (handle from p3).
    p2: [f64; DIMS_MAX],
    /// End point (on the curve).
    p3: [f64; DIMS_MAX],
}

// ============================================================================
// Bernstein Polynomials
// ============================================================================

/// Compute all four Bernstein polynomial values with shared intermediate calculations.
/// Avoids redundant computation of (1-u), u^2, etc. when all four values are needed.
#[inline]
fn bernstein_all(u: f64) -> (f64, f64, f64, f64) {
    let s = 1.0 - u;
    let ss = s * s;
    let uu = u * u;
    let us3 = 3.0 * u * s;

    let b1 = us3 * s;                      // 3 * u * (1-u)^2
    let b2 = us3 * u;                      // 3 * u^2 * (1-u)
    let b0_plus_b1 = ss * (1.0 + 2.0 * u);
    let b2_plus_b3 = uu * (3.0 - 2.0 * u);

    (b1, b2, b0_plus_b1, b2_plus_b3)
}

// ============================================================================
// Circular Fallback Helpers
// ============================================================================

/// Return a scale value based on circular arc approximation.
///
/// This works by placing each end-point on an imaginary circle,
/// the placement on the circle is based on the tangent vectors,
/// where larger differences in tangent angle cover a larger part of the circle.
///
/// Returns the scale representing how much larger the distance around the circle is.
fn points_calc_circumference_factor(
    tan_l: &[f64],
    tan_r: &[f64],
) -> f64 {
    use crate::math_vector::{len_vnvn, len_negated_vnvn};
    use std::f64;

    let dot = dot_vnvn(tan_l, tan_r);

    let len_tangent = if dot < 0.0 { len_vnvn } else { len_negated_vnvn } (tan_l, tan_r);
    if len_tangent > f64::EPSILON {
        // Only clamp to avoid precision error.
        let angle = ((-dot.abs()).max(-1.0)).acos();
        // Angle may be less than the length when the
        // tangents define >180 degrees of the circle,
        // (tangents that point away from each other).
        //
        // We could try support this but will likely cause
        // extreme >1 scales which could cause other issues.

        // assert(angle >= len_tangent);
        let factor = angle / len_tangent;
        debug_assert!(factor < (f64::consts::PI / 2.0) + (f64::EPSILON * 10.0));
        factor
    } else {
        // Tangents are exactly aligned (think two opposite sides of a circle).
        f64::consts::PI / 2.0
    }
}

/// Return the handle scale factor for points on a perfect circle.
///
/// Note: the return value will need to be multiplied by 1.3... for correct results.
fn points_calc_circle_tangent_factor(
    tan_l: &[f64],
    tan_r: &[f64],
) -> Option<f64> {
    let eps = 1e-8;
    let tan_dot = dot_vnvn(tan_l, tan_r);
    if tan_dot > 1.0 - eps {
        // No angle difference (use fallback, length won't make any difference).
        Some((1.0 / 3.0) * 0.75)
    } else if tan_dot < -1.0 + eps {
        // Parallel tangents (half-circle).
        Some(1.0 / 2.0)
    } else {
        // Non-aligned tangents, calculate handle length.
        let angle = tan_dot.acos() / 2.0;

        // Could also use 'angle_sin = len_vnvn(tan_l, tan_r) / 2.0'.
        let angle_sin = angle.sin();
        let angle_cos = angle.cos();
        Some(((1.0 - angle_cos) / (angle_sin * 2.0)) / angle_sin)
    }
}

/// Calculate the handle scale using circular arc approximation.
///
/// Serves as a best-guess fallback when the least-squares solution fails.
fn points_calc_cubic_scale(
    v_l: &[f64],
    v_r: &[f64],
    tan_l: &[f64],
    tan_r: &[f64],
    coords_length: f64,
) -> Option<f64> {
    if let Some(len_circle_factor) = points_calc_circle_tangent_factor(tan_l, tan_r) {
        let len_direct = len_vnvn(v_l, v_r);

        // If this curve is a circle, this value doesn't need modification.
        let len_circle_handle = len_direct * (len_circle_factor / 0.75);

        // Scale by the difference from the circumference distance.
        let len_circle = len_direct * points_calc_circumference_factor(tan_l, tan_r);
        let mut scale_handle = coords_length / len_circle;

        // Could investigate an accurate calculation here,
        // though this gives close results.
        scale_handle = ((scale_handle - 1.0) * 1.75) + 1.0;

        scale_handle *= len_circle_handle;

        if scale_handle.is_finite() {
            return Some(scale_handle);
        }
    }
    None
}

// ============================================================================
// Cubic Solvers
// ============================================================================

/// Simple fallback: calculate handles based on 1/3 of endpoint distance.
///
/// This is used as a baseline when more sophisticated methods fail.
fn cubic_from_points_fallback(
    points: &[f64],
    dims: usize,
    tan_l: &[f64],
    tan_r: &[f64],
) -> Cubic {
    use crate::math_vector::{madd_vnvn_fl, msub_vnvn_fl};

    let points_len = points.len() / dims;
    let p0 = &points[0..dims];
    let p3 = &points[(points_len - 1) * dims..];
    let alpha = len_vnvn(p0, p3) / 3.0;

    let mut p0_arr = [0.0; DIMS_MAX];
    let mut p3_arr = [0.0; DIMS_MAX];
    p0_arr[..dims].copy_from_slice(p0);
    p3_arr[..dims].copy_from_slice(p3);

    Cubic {
        p0: p0_arr,
        p1: msub_vnvn_fl(p0, tan_l, alpha, dims),
        p2: madd_vnvn_fl(p3, tan_r, alpha, dims),
        p3: p3_arr,
    }
}

/// Use least-squares method to find Bezier control points for region.
///
/// Returns (cubic, use_clamp) where use_clamp indicates if handle clamping should be applied.
/// Returns None for the cubic if the solution produces invalid (negative) alpha values.
fn cubic_from_points(
    points: &[f64],
    dims: usize,
    points_coords_length: f64,
    u_prime: &[f64],
    tan_l: &[f64],
    tan_r: &[f64],
) -> (Cubic, bool) {
    use crate::math_vector::{mul_vn_fl, madd_vnvn_fl, msub_vnvn_fl, is_almost_zero};

    let points_len = points.len() / dims;
    let p0 = &points[0..dims];
    let p3 = &points[(points_len - 1) * dims..];

    let (alpha_l, alpha_r) = {
        let mut x: [f64; 2] = [0.0, 0.0];
        let mut c: [[f64; 2]; 2] = [[0.0, 0.0], [0.0, 0.0]];

        for (i, u) in u_prime.iter().enumerate() {
            let pt = &points[i * dims..(i + 1) * dims];
            let (b1, b2, b0_plus_b1, b2_plus_b3) = bernstein_all(*u);
            let a0 = mul_vn_fl(tan_l, b1, dims);
            let a1 = mul_vn_fl(tan_r, b2, dims);

            // Inline dot product.
            for j in 0..dims {
                let tmp = (pt[j] - (p0[j] * b0_plus_b1)) + (p3[j] * b2_plus_b3);

                x[0] += a0[j] * tmp;
                x[1] += a1[j] * tmp;

                c[0][0] += a0[j] * a0[j];
                c[0][1] += a0[j] * a1[j];
                c[1][1] += a1[j] * a1[j];
            }

            c[1][0] = c[0][1];
        }

        let det_c0_c1 = {
            let tmp = c[0][0] * c[1][1] - c[0][1] * c[1][0];
            if !is_almost_zero(tmp) {
                tmp
            } else {
                c[0][0] * c[1][1] * 10e-12
            }
        };
        let det_c_0x = x[1] * c[0][0] - x[0] * c[0][1];
        let det_x_c1 = x[0] * c[1][1] - x[1] * c[0][1];

        let alpha_l = det_x_c1 / det_c0_c1;
        let alpha_r = det_c_0x / det_c0_c1;

        // May still divide-by-zero, check below will catch NaN values.
        (alpha_l, alpha_r)
    };

    let mut use_clamp = true;

    // Flip check to catch NaN values.
    let (alpha_l, alpha_r) = if !(alpha_l >= 0.0) || !(alpha_r >= 0.0) {
        // Use circular fallback.
        let alpha_test = if let Some(scale) = points_calc_cubic_scale(p0, p3, tan_l, tan_r, points_coords_length) {
            scale
        } else {
            len_vnvn(p0, p3) / 3.0
        };
        // Skip clamping when we're using default handles.
        use_clamp = false;
        (alpha_test, alpha_test)
    } else {
        (alpha_l, alpha_r)
    };

    let mut p0_arr = [0.0; DIMS_MAX];
    let mut p3_arr = [0.0; DIMS_MAX];
    p0_arr[..dims].copy_from_slice(p0);
    p3_arr[..dims].copy_from_slice(p3);

    let cubic = Cubic {
        p0: p0_arr,
        p1: msub_vnvn_fl(p0, tan_l, alpha_l, dims),
        p2: madd_vnvn_fl(p3, tan_r, alpha_r, dims),
        p3: p3_arr,
    };

    (cubic, use_clamp)
}

/// Offset-based fallback: use perpendicular distance from the line-segment.
///
/// Uses the maximum perpendicular distance of points from the line
/// between endpoints to determine handle lengths. Can do a 'perfect'
/// reversal of subdivision when the curve has symmetrical handles.
fn cubic_from_points_offset_fallback(
    points: &[f64],
    dims: usize,
    tan_l: &[f64],
    tan_r: &[f64],
) -> Cubic {
    use crate::math_vector::{
        madd_vnvn_fl, msub_vnvn_fl,
        negated_vn,
        normalized_vnvn,
        normalized_vn,
        project_plane_vnvn_normalized,
        project_vnvn_normalized,
    };
    use std::f64;

    let points_len = points.len() / dims;
    let p0 = &points[0..dims];
    let p3 = &points[(points_len - 1) * dims..];

    let dir_dist = len_vnvn(p0, p3);
    let dir_unit = normalized_vnvn(p3, p0, dims);
    // Note that normalizing output here is only for better accuracy, not essential.
    let a0 = normalized_vn(&project_plane_vnvn_normalized(tan_l, &dir_unit[..dims], dims), dims);
    let a1 = negated_vn(&normalized_vn(&project_plane_vnvn_normalized(tan_r, &dir_unit[..dims], dims), dims), dims);

    let mut dists: [f64; 2] = [0.0, 0.0];

    for i in 1..(points_len - 1) {
        let pt = &points[i * dims..(i + 1) * dims];
        let sub0 = sub_vnvn(p0, pt, dims);
        let tmp0 = project_vnvn_normalized(&sub0[..dims], &a0[..dims], dims);
        dists[0] = dists[0].max(dot_vnvn(&tmp0[..dims], &a0[..dims]));
        let sub1 = sub_vnvn(p0, pt, dims);
        let tmp1 = project_vnvn_normalized(&sub1[..dims], &a1[..dims], dims);
        dists[1] = dists[1].max(dot_vnvn(&tmp1[..dims], &a1[..dims]));
    }

    // The value of 'dists[..] / 0.75' is the length to use when the tangents
    // are perpendicular to the direction defined by the two points.
    //
    // Project tangents onto these perpendicular lengths.
    // Note that this can cause divide by zero in the case of collinear tangents.
    // The limits check afterwards accounts for this.
    let div_l = dot_vnvn(&tan_l[..dims], &a0[..dims]).abs();
    let div_r = dot_vnvn(&tan_r[..dims], &a1[..dims]).abs();

    let mut alpha_l = if div_l > 0.0 { (dists[0] / 0.75) / div_l } else { f64::INFINITY };
    let mut alpha_r = if div_r > 0.0 { (dists[1] / 0.75) / div_r } else { f64::INFINITY };

    if !(alpha_l > 0.0) || (alpha_l > dists[0] + dir_dist) {
        alpha_l = dir_dist / 3.0;
    }
    if !(alpha_r > 0.0) || (alpha_r > dists[1] + dir_dist) {
        alpha_r = dir_dist / 3.0;
    }

    let mut p0_arr = [0.0; DIMS_MAX];
    let mut p3_arr = [0.0; DIMS_MAX];
    p0_arr[..dims].copy_from_slice(p0);
    p3_arr[..dims].copy_from_slice(p3);

    Cubic {
        p0: p0_arr,
        p1: msub_vnvn_fl(p0, tan_l, alpha_l, dims),
        p2: madd_vnvn_fl(p3, tan_r, alpha_r, dims),
        p3: p3_arr,
    }
}


/// Use Newton-Raphson iteration to find better root.
///
/// * `cubic` - Current fitted curve.
/// * `p` - Point to test against.
/// * `u` - Parameter value for `p`.
/// * `dims` - Number of dimensions.
///
/// Note: return value may be `nan` caller must check for this.
fn cubic_find_root(
    cubic: &Cubic,
    p: &[f64],
    u: f64,
    dims: usize,
) -> f64 {
    // Newton-Raphson Method.
    let (point, q1_u, q2_u) = cubic_calc_point_speed_accel(cubic, u, dims);
    let q0_u = sub_vnvn(&point[..dims], p, dims);

    // May divide-by-zero, caller must check for that case.
    // u - (q0_u * q1_u) / (q1_u.length_squared() + q0_u * q2_u)
    u - dot_vnvn(&q0_u[..dims], &q1_u[..dims]) / (dot_vnvn(&q1_u[..dims], &q1_u[..dims]) + dot_vnvn(&q0_u[..dims], &q2_u[..dims]))
}

/// Given set of points and their parameterization, try to find a better parameterization.
///
/// Uses Newton-Raphson iteration on each point to find parameter values that
/// minimize distance to the curve. Returns false if the reparameterization
/// produces invalid values (non-finite, out of [0,1] range, or unsorted).
fn cubic_reparameterize(
    cubic: &Cubic,
    points: &[f64],
    dims: usize,
    u_prime_src: &[f64],
    u_prime_dst: &mut [f64]
) -> bool {
    let points_len = points.len() / dims;
    debug_assert!(points_len == u_prime_src.len());
    debug_assert!(points_len == u_prime_dst.len());


    // Recalculate the values of u[] based on the Newton-Raphson method.
    for i in 0..points_len {
        let pt = &points[i * dims..(i + 1) * dims];
        u_prime_dst[i] = cubic_find_root(cubic, pt, u_prime_src[i], dims);
        if !u_prime_dst[i].is_finite() {
            return false;
        }
    }

    // we can safely unwrap here because nan/inf's are caught above
    u_prime_dst.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if (u_prime_dst[0] < 0.0) ||
       (u_prime_dst[points_len - 1] > 1.0)
    {
        return false;
    }

    debug_assert!(u_prime_dst[0] >= 0.0);
    debug_assert!(u_prime_dst[u_prime_dst.len() - 1] <= 1.0);
    true
}

/// Calculate normalized arc-length parameterization for points.
///
/// Returns a tuple of (u, total_length) where:
/// - `u` is a vector of parameter values in [0, 1] for each point
/// - `total_length` is the total arc length of the polyline
fn points_calc_coord_length(
    points: &[f64],
    dims: usize,
    points_length_cache: &[f64],
) -> (Vec<f64>, f64) {
    let points_len = points.len() / dims;
    let mut u: Vec<f64> = Vec::with_capacity(points_len);
    u.push(0.0);

    let mut l_prev = 0.0;
    for i in 1..points_len {
        let pt_prev = &points[(i - 1) * dims..i * dims];
        let pt = &points[i * dims..(i + 1) * dims];
        let l = points_length_cache[i];
        debug_assert!((len_vnvn(pt, pt_prev) - l).abs() < 1e-10);
        let l_curr = l + l_prev;
        u.push(l_curr);
        l_prev = l_curr;
    }

    debug_assert!(u.len() == points_len);

    let w = u[u.len() - 1];
    let w_inv = 1.0 / w;
    for u_step in &mut u[1..] {
        *u_step *= w_inv;
    }

    (u, w)
}

/// Evaluate the cubic bezier curve at parameter t using de Casteljau's algorithm.
fn cubic_calc_point(
    cubic: &Cubic, t: f64, dims: usize,
) -> [f64; DIMS_MAX] {
    let p0 = &cubic.p0;
    let p1 = &cubic.p1;
    let p2 = &cubic.p2;
    let p3 = &cubic.p3;
    let s = 1.0 - t;
    let mut v_out = [0.0; DIMS_MAX];
    for j in 0..dims {
        let p01 = (p0[j] * s) + (p1[j] * t);
        let p12 = (p1[j] * s) + (p2[j] * t);
        let p23 = (p2[j] * s) + (p3[j] * t);
        v_out[j] = ((((p01 * s) + (p12 * t))) * s) +
                   ((((p12 * s) + (p23 * t))) * t);
    }
    v_out
}

/// Compute point, first derivative (speed), and second derivative (acceleration) in one pass.
/// Combines #cubic_calc_point, #cubic_calc_speed, and #cubic_calc_acceleration to share
/// intermediate values and reduce redundant control point access.
#[inline]
fn cubic_calc_point_speed_accel(
    cubic: &Cubic, t: f64, dims: usize,
) -> ([f64; DIMS_MAX], [f64; DIMS_MAX], [f64; DIMS_MAX]) {
    let p0 = &cubic.p0;
    let p1 = &cubic.p1;
    let p2 = &cubic.p2;
    let p3 = &cubic.p3;
    let s = 1.0 - t;
    let ss = s * s;
    let tt = t * t;
    let st2 = 2.0 * s * t;

    let mut r_point = [0.0; DIMS_MAX];
    let mut r_speed = [0.0; DIMS_MAX];
    let mut r_accel = [0.0; DIMS_MAX];

    for j in 0..dims {
        // Control point differences, computed once.
        let d01 = p1[j] - p0[j];
        let d12 = p2[j] - p1[j];
        let d23 = p3[j] - p2[j];

        // Point via de Casteljau's algorithm.
        let p01 = p0[j] + d01 * t;
        let p12 = p1[j] + d12 * t;
        let p23 = p2[j] + d23 * t;
        let p012 = p01 * s + p12 * t;
        let p123 = p12 * s + p23 * t;
        r_point[j] = p012 * s + p123 * t;

        // First derivative: 3 * ((d01)*s^2 + 2*(d12)*s*t + (d23)*t^2).
        r_speed[j] = 3.0 * (d01 * ss + d12 * st2 + d23 * tt);

        // Second derivative: 6 * ((d12 - d01)*s + (d23 - d12)*t).
        r_accel[j] = 6.0 * ((d12 - d01) * s + (d23 - d12) * t);
    }

    (r_point, r_speed, r_accel)
}

/// Error metrics from fitting a cubic to a set of points.
#[derive(Clone, Copy)]
struct FitError {
    /// Maximum squared distance from any point to the fitted curve.
    pub max_sq: f64,
    /// Index of the point with maximum error (potential split point).
    pub index: usize,
}

/// Returns a 'measure' of the maximum distance (squared) of the points
/// from the corresponding cubic(u[]) points.
///
/// Returns the maximum squared error and the index of the point with that error.
/// The index can be used as a split point if the error exceeds the threshold.
fn cubic_calc_error(
    cubic: &Cubic,
    points: &[f64],
    dims: usize,
    u: &[f64],
) -> FitError {
    let points_len = points.len() / dims;
    let mut error_max_sq = -1.0;

    // No need to measure first & last points (they are on the curve by construction).
    let mut index = 1;
    let mut error_index = 1;
    for i in 1..(points_len - 1) {
        let pt_real = &points[i * dims..(i + 1) * dims];
        let u_step = u[i];
        let pt_eval = cubic_calc_point(cubic, u_step, dims);
        let err_sq = len_squared_vnvn(pt_real, &pt_eval[..dims]);
        // Use >= to match C behavior: pick the last point with max error.
        if err_sq >= error_max_sq {
            error_max_sq = err_sq;
            error_index = index;
        }
        index += 1;
    }

    debug_assert!(error_max_sq != -1.0);
    FitError {
        max_sq: error_max_sq,
        index: error_index,
    }
}

/// Like #cubic_calc_error but return None
/// in the case we can't improve on `error_max_sq_limit`.
fn cubic_calc_error_limit(
    cubic: &Cubic,
    points: &[f64],
    dims: usize,
    u: &[f64],
    error_max_sq_limit: f64,
) -> Option<FitError> {
    let points_len = points.len() / dims;
    let mut error_max_sq = -1.0;

    // no need to measure first & last points
    let mut index = 1;
    let mut error_index = 1;
    for i in 1..(points_len - 1) {
        let pt_real = &points[i * dims..(i + 1) * dims];
        let u_step = u[i];
        let pt_eval = cubic_calc_point(cubic, u_step, dims);
        let err_sq = len_squared_vnvn(pt_real, &pt_eval[..dims]);
        // Use >= to match C behavior.
        if err_sq >= error_max_sq_limit {
            return None;
        } else if err_sq >= error_max_sq {
            error_max_sq = err_sq;
            error_index = index;
        }
        index += 1;
    }

    debug_assert!(error_max_sq != -1.0);
    Some(FitError {
        max_sq: error_max_sq,
        index: error_index,
    })
}

/// Calculate a weighted center that compensates for non-uniform point spacing.
///
/// Each point is weighted by the sum of distances to its neighbors,
/// giving more influence to points in sparse regions. This provides
/// a more representative center for handle clamping calculations.
fn points_calc_center_weighted(
    points: &[f64],
    dims: usize,
) -> [f64; DIMS_MAX] {
    let points_len = points.len() / dims;
    let mut center = [0.0; DIMS_MAX];
    let mut w_tot = 0.0;

    let pt_prev_start = &points[(points_len - 2) * dims..(points_len - 1) * dims];
    let pt_curr_start = &points[(points_len - 1) * dims..];

    let mut w_prev = len_vnvn(pt_prev_start, pt_curr_start);
    let mut pt_curr = pt_curr_start;

    for i_next in 0..points_len {
        let pt_next = &points[i_next * dims..(i_next + 1) * dims];
        let w_next = len_vnvn(pt_curr, pt_next);
        let w = w_prev + w_next;
        w_tot += w;

        for j in 0..dims {
            center[j] += pt_curr[j] * w;
        }

        w_prev = w_next;
        pt_curr = pt_next;
    }

    if w_tot != 0.0 {
        let w_inv = 1.0 / w_tot;
        for j in 0..dims {
            center[j] *= w_inv;
        }
    }

    center
}

/// Apply handle clamping to prevent extreme handle values.
/// Clamps handles to be within 3x the maximum distance from any point to the weighted center.
fn cubic_apply_handle_clamping(
    cubic: &mut Cubic,
    points: &[f64],
    dims: usize,
    tan_l: &[f64],
    tan_r: &[f64],
    points_length: f64,
) {
    let points_len = points.len() / dims;
    let center = points_calc_center_weighted(points, dims);

    /// Maximum handle distance as a multiple of the point cloud radius.
    const CLAMP_SCALE: f64 = 3.0;

    // Find max distance squared from center to any point, scaled by clamp_scale.
    let mut dist_sq_max: f64 = 0.0;
    for i in 0..points_len {
        let pt = &points[i * dims..(i + 1) * dims];
        let mut dist_sq_test: f64 = 0.0;
        for j in 0..dims {
            dist_sq_test += sq((pt[j] - center[j]) * CLAMP_SCALE);
        }
        dist_sq_max = dist_sq_max.max(dist_sq_test);
    }

    let mut p1_dist_sq = len_squared_vnvn(&center[..dims], &cubic.p1[..dims]);
    let mut p2_dist_sq = len_squared_vnvn(&center[..dims], &cubic.p2[..dims]);

    // If either handle exceeds the limit, fall back to a simpler calculation.
    if p1_dist_sq > dist_sq_max || p2_dist_sq > dist_sq_max {
        // Try circular fallback scale.
        let alpha_test = if let Some(scale) = points_calc_cubic_scale(
            &cubic.p0[..dims], &cubic.p3[..dims], tan_l, tan_r, points_length)
        {
            scale
        } else {
            len_vnvn(&cubic.p0[..dims], &cubic.p3[..dims]) / 3.0
        };

        // Recalculate handles with fallback alpha.
        for j in 0..dims {
            cubic.p1[j] = cubic.p0[j] - tan_l[j] * alpha_test;
            cubic.p2[j] = cubic.p3[j] + tan_r[j] * alpha_test;
        }

        p1_dist_sq = len_squared_vnvn(&center[..dims], &cubic.p1[..dims]);
        p2_dist_sq = len_squared_vnvn(&center[..dims], &cubic.p2[..dims]);
    }

    // Clamp handles within the 3x radius.
    if p1_dist_sq > dist_sq_max {
        let scale = (dist_sq_max.sqrt()) / (p1_dist_sq.sqrt());
        for j in 0..dims {
            cubic.p1[j] = center[j] + (cubic.p1[j] - center[j]) * scale;
        }
    }
    if p2_dist_sq > dist_sq_max {
        let scale = (dist_sq_max.sqrt()) / (p2_dist_sq.sqrt());
        for j in 0..dims {
            cubic.p2[j] = center[j] + (cubic.p2[j] - center[j]) * scale;
        }
    }
}

/// Attempt to fit a cubic bezier curve to a set of points.
///
/// This function matches the C algorithm exactly:
/// 1. Compute parameterization and try least-squares (with circular fallback if LS fails).
/// 2. If error < threshold: return early.
/// 3. Try simple fallback (1/3 distance) - only if above threshold.
/// 4. If error < threshold: return early.
/// 5. Try offset fallback - only if above threshold.
/// 6. If error < threshold: return early.
/// 7. Iteration loop (max 4) - only if above threshold, with early exit on success.
///
/// Returns (cubic, error, threshold_met).
fn fit_cubic_to_points(
    points: &[f64],
    dims: usize,
    points_length_cache: &[f64],
    tan_l: &[f64],
    tan_r: &[f64],
    error_threshold_sq: f64,
) -> (Cubic, FitError, bool) {
    // Maximum Newton-Raphson iterations for parameter refinement.
    let iteration_max = 4;
    let points_len = points.len() / dims;

    // Special case: 2 points - create a simple cubic with handles at 1/3 distance.
    if points_len == 2 {
        let p0 = &points[0..dims];
        let p3 = &points[dims..];
        let dist = len_vnvn(p0, p3) / 3.0;

        let mut cubic = Cubic {
            p0: [0.0; DIMS_MAX],
            p1: [0.0; DIMS_MAX],
            p2: [0.0; DIMS_MAX],
            p3: [0.0; DIMS_MAX],
        };
        cubic.p0[..dims].copy_from_slice(p0);
        cubic.p3[..dims].copy_from_slice(p3);
        for j in 0..dims {
            cubic.p1[j] = p0[j] - tan_l[j] * dist;
            cubic.p2[j] = p3[j] + tan_r[j] * dist;
        }

        // Error is 0 for a 2-point curve (perfect fit).
        let error = FitError { max_sq: 0.0, index: 0 };
        return (cubic, error, true);
    }

    assert!(points_len > 2);

    // Step 1: Compute parameterization.
    let (mut u, points_length) = points_calc_coord_length(points, dims, points_length_cache);

    // Step 2: Call least-squares (with circular fallback if LS fails), apply clamping.
    // This matches C's cubic_from_points() behavior.
    let (mut cubic_best, use_clamp) = cubic_from_points(
        points, dims, points_length, &u, tan_l, tan_r);
    if use_clamp {
        cubic_apply_handle_clamping(&mut cubic_best, points, dims, tan_l, tan_r, points_length);
    }

    let mut error_best = cubic_calc_error(&cubic_best, points, dims, &u);

    // Early exit if already within threshold.
    if error_best.max_sq < error_threshold_sq {
        return (cubic_best, error_best, true);
    }

    // Step 3: Try simple fallback (cubic_from_points_fallback in C).
    // C uses simple 1/3 fallback here, NOT the full cubic_from_points.
    let split_index;
    {
        let cubic_fallback = cubic_from_points_fallback(points, dims, tan_l, tan_r);
        let error_fallback = cubic_calc_error(&cubic_fallback, points, dims, &u);

        // C: "Intentionally use the newly calculated 'split_index',
        // even if the 'error_max_sq_test' is worse."
        split_index = error_fallback.index;

        if error_best.max_sq > error_fallback.max_sq {
            error_best = error_fallback;
            cubic_best = cubic_fallback;
        }
    }

    if error_best.max_sq < error_threshold_sq {
        error_best.index = split_index;
        return (cubic_best, error_best, true);
    }

    // Step 4: Try offset fallback.
    {
        let cubic_offset = cubic_from_points_offset_fallback(points, dims, tan_l, tan_r);
        // C uses #cubic_calc_error_simple which returns early if error >= limit.
        if let Some(error_offset) = cubic_calc_error_limit(
            &cubic_offset, points, dims, &u, error_best.max_sq)
        {
            error_best = error_offset;
            cubic_best = cubic_offset;
        }
    }

    // Update split_index in error_best.
    error_best.index = split_index;

    if error_best.max_sq < error_threshold_sq {
        return (cubic_best, error_best, true);
    }

    // Step 5: Iteration loop (only if still above threshold).
    let mut cubic_for_reparam = cubic_best.clone();
    let mut u_prime: Vec<f64> = vec![0.0; u.len()];

    for _iter in 0..iteration_max {
        if !cubic_reparameterize(&cubic_for_reparam, points, dims, &u, &mut u_prime) {
            break;
        }

        // Call least-squares (with circular fallback if LS fails), apply clamping.
        // This matches C's cubic_from_points() behavior.
        let (mut cubic_test, use_clamp) = cubic_from_points(
            points, dims, points_length, &u_prime, tan_l, tan_r);
        if use_clamp {
            cubic_apply_handle_clamping(&mut cubic_test, points, dims, tan_l, tan_r, points_length);
        }
        let error_test = cubic_calc_error(&cubic_test, points, dims, &u_prime);

        // Always update cubic_for_reparam for next iteration (matching C behavior).
        cubic_for_reparam = cubic_test.clone();

        if error_best.max_sq > error_test.max_sq {
            error_best = error_test;
            cubic_best = cubic_test;
        }

        // Early exit if below threshold.
        if error_best.max_sq < error_threshold_sq {
            return (cubic_best, error_best, true);
        }

        std::mem::swap(&mut u, &mut u_prime);
    }

    // Threshold not met.
    (cubic_best, error_best, false)
}

/// Fit a cubic bezier curve to a set of points.
///
/// Returns ((error_sq, split_index), handle_left, handle_right, threshold_met) where:
/// - `error_sq` is the maximum squared error from any point to the curve
/// - `split_index` is the index of the point with maximum error
/// - `handle_left` and `handle_right` are the bezier control points (p1, p2)
/// - `threshold_met` is true if error_sq < error_threshold_sq
pub fn curve_fit_cubic_to_points_single(
    points: &[f64],
    dims: usize,
    points_length_cache: &[f64],
    tan_l: &[f64],
    tan_r: &[f64],
    error_threshold_sq: f64,
) -> ((f64, usize), [f64; DIMS_MAX], [f64; DIMS_MAX], bool) {
    let (cubic, fit_error, threshold_met) = fit_cubic_to_points(
        points,
        dims,
        points_length_cache,
        tan_l, tan_r,
        error_threshold_sq);

    ((fit_error.max_sq, fit_error.index), cubic.p1, cubic.p2, threshold_met)
}

// =============================================================================
// Split Point Calculation Functions
// =============================================================================
//
// These functions find candidate split points for curve subdivision.
// They operate on raw point arrays with index ranges, matching the C API.

/// Invalid split point marker.
pub const SPLIT_POINT_INVALID: usize = usize::MAX;

/// Find the split point with maximum perpendicular distance from the line-segment.
///
/// # Returns
/// Index of the split point, or SPLIT_POINT_INVALID if none found.
pub fn split_point_find_max_distance(
    points: &[f64],
    points_len: usize,
    index_l: usize,
    index_r: usize,
    dims: usize,
) -> usize {
    use crate::math_vector::{
        normalize_vn, project_plane_vnvn_normalized, len_squared_vn,
    };

    let mut split_point = SPLIT_POINT_INVALID;
    let mut split_point_dist_best: f64 = -f64::MAX;

    let offset = &points[index_l * dims..(index_l + 1) * dims];

    // Direction along the segment (line from `index_l` to `index_r`).
    let mut v_segment = [0.0; DIMS_MAX];
    for j in 0..dims {
        v_segment[j] = points[index_l * dims + j] - points[index_r * dims + j];
    }
    normalize_vn(&mut v_segment[..dims]);

    // Iterate from `index_l + 1` to `index_r - 1` (exclusive of endpoints).
    let mut i_curr = index_l;
    loop {
        // Advance to next index, wrapping at points_len.
        i_curr = (i_curr + 1) % points_len;

        if i_curr == index_r {
            break;
        }

        let mut v_offset = [0.0; DIMS_MAX];
        for j in 0..dims {
            v_offset[j] = points[i_curr * dims + j] - offset[j];
        }
        let v_proj = project_plane_vnvn_normalized(&v_offset[..dims], &v_segment[..dims], dims);

        let dist_sq = len_squared_vn(&v_proj[..dims]);
        if dist_sq > split_point_dist_best {
            split_point_dist_best = dist_sq;
            split_point = i_curr;
        }
    }

    split_point
}

/// Find the split point based on sign change of perpendicular distance.
///
/// This finds where the curve crosses the line between the two endpoints,
/// selecting the crossing point with the largest perpendicular distance.
///
/// # Returns
/// Index of the split point, or SPLIT_POINT_INVALID if none found.
///
/// NOTE: This operation is symmetrical (reversing the curve produces the same split point).
pub fn split_point_find_sign_change(
    points: &[f64],
    points_len: usize,
    index_l: usize,
    index_r: usize,
    dims: usize,
) -> usize {
    use crate::math_vector::{
        normalize_vn, project_plane_vnvn_normalized, len_squared_vn,
    };

    let mut split_point = SPLIT_POINT_INVALID;
    let mut split_point_dist_best: f64 = -f64::MAX;

    let offset = &points[index_l * dims..(index_l + 1) * dims];

    // Direction along the line.
    let mut v_line = [0.0; DIMS_MAX];
    for j in 0..dims {
        v_line[j] = points[index_l * dims + j] - points[index_r * dims + j];
    }
    normalize_vn(&mut v_line[..dims]);

    // Reference perpendicular direction.
    let mut v_ref = [0.0; DIMS_MAX];
    let mut have_reference = false;

    let mut i_curr = index_l;
    let mut i_best = index_l;
    let mut best_signed_dist: f64 = 0.0;
    let mut best_dist_sq: f64 = 0.0;

    loop {
        i_curr = (i_curr + 1) % points_len;

        if i_curr == index_r {
            break;
        }

        let mut v_offset = [0.0; DIMS_MAX];
        for j in 0..dims {
            v_offset[j] = points[i_curr * dims + j] - offset[j];
        }
        let v_proj = project_plane_vnvn_normalized(&v_offset[..dims], &v_line[..dims], dims);

        let dist_sq = len_squared_vn(&v_proj[..dims]);

        // Establish reference direction from first point with significant deviation.
        if !have_reference && dist_sq > 1e-12 {
            let len_inv = 1.0 / dist_sq.sqrt();
            for j in 0..dims {
                v_ref[j] = v_proj[j] * len_inv;
            }
            have_reference = true;
        }

        // Signed distance is the dot product with the reference direction.
        let signed_dist = if have_reference {
            dot_vnvn(&v_proj[..dims], &v_ref[..dims])
        } else {
            0.0
        };

        // Check for sign change (curve crossing the line).
        if best_signed_dist * signed_dist < 0.0 {
            if dist_sq > split_point_dist_best {
                split_point_dist_best = dist_sq;
                split_point = i_curr;
            }
            if best_dist_sq > split_point_dist_best {
                split_point_dist_best = best_dist_sq;
                split_point = i_best;
            }
        }

        best_signed_dist = signed_dist;
        best_dist_sq = dist_sq;
        i_best = i_curr;
    }

    split_point
}

/// Find the split point with maximum projection onto a given axis.
///
/// Used for corner detection - finds the point that extends furthest
/// in the direction the corner "points".
///
/// # Returns
/// Index of the split point, or SPLIT_POINT_INVALID if none found.
pub fn split_point_find_max_on_axis(
    points: &[f64],
    points_len: usize,
    index_l: usize,
    index_r: usize,
    axis: &[f64],
    dims: usize,
) -> usize {
    let mut split_point = SPLIT_POINT_INVALID;
    let mut split_point_proj_best: f64 = -f64::MAX;

    // Iterate from index_l+1 to index_r-1 (exclusive of endpoints).
    let mut i_curr = index_l;
    loop {
        // Advance to next index, wrapping at points_len.
        i_curr = (i_curr + 1) % points_len;

        if i_curr == index_r {
            break;
        }

        let proj = dot_vnvn(axis, &points[i_curr * dims..(i_curr + 1) * dims]);
        if proj > split_point_proj_best {
            split_point_proj_best = proj;
            split_point = i_curr;
        }
    }

    split_point
}

/// Find the split point based on inflection (curvature sign change).
///
/// This finds where the curve changes from curving one way to curving the other,
/// selecting the inflection point with the largest curvature magnitude.
///
/// # Returns
/// Index of the split point, or SPLIT_POINT_INVALID if none found.
///
/// NOTE: This operation is symmetrical (reversing the curve produces the same split point).
pub fn split_point_find_inflection(
    points: &[f64],
    points_len: usize,
    index_l: usize,
    index_r: usize,
    dims: usize,
) -> usize {
    use crate::math_vector::len_squared_vn;

    let mut split_point = SPLIT_POINT_INVALID;
    let mut split_point_accel_best: f64 = -f64::MAX;

    // Calculate span length to check if we have enough points.
    let span_len = if index_l <= index_r {
        index_r - index_l + 1
    } else {
        index_r + points_len - index_l + 1
    };

    // Need at least 4 points to detect inflection.
    if span_len < 4 {
        return split_point;
    }

    let mut v_accel = [0.0; DIMS_MAX];
    let mut v_ref = [0.0; DIMS_MAX];

    let mut i_prev = index_l;
    let mut i_curr = (index_l + 1) % points_len;

    let mut i_best = i_curr;
    let mut best_signed_accel: f64 = 0.0;
    let mut best_accel_sq: f64 = 0.0;
    let mut have_reference = false;

    // Safety limit to prevent infinite loops.
    let mut iter_limit = points_len + 1;

    while i_curr != index_r {
        iter_limit -= 1;
        if iter_limit == 0 {
            break;
        }

        let i_next = (i_curr + 1) % points_len;

        // Compute acceleration: accel = p[i+1] - 2*p[i] + p[i-1]
        // This is the discrete second derivative (finite difference) on raw input points.
        // Unlike #cubic_calc_acceleration which computes the parametric second derivative
        // of a Bezier curve at parameter t, this operates directly on sampled points.
        for j in 0..dims {
            v_accel[j] = points[i_next * dims + j]
                       - 2.0 * points[i_curr * dims + j]
                       + points[i_prev * dims + j];
        }

        let accel_sq = len_squared_vn(&v_accel[..dims]);

        // Establish reference direction.
        if !have_reference && accel_sq > 1e-12 {
            let len_inv = 1.0 / accel_sq.sqrt();
            for j in 0..dims {
                v_ref[j] = v_accel[j] * len_inv;
            }
            have_reference = true;
        }

        let signed_accel = if have_reference {
            dot_vnvn(&v_accel[..dims], &v_ref[..dims])
        } else {
            0.0
        };

        // Check for sign change (inflection point).
        if best_signed_accel * signed_accel < 0.0 {
            if accel_sq > split_point_accel_best {
                split_point_accel_best = accel_sq;
                split_point = i_curr;
            }
            if best_accel_sq > split_point_accel_best {
                split_point_accel_best = best_accel_sq;
                split_point = i_best;
            }
        }

        best_signed_accel = signed_accel;
        best_accel_sq = accel_sq;
        i_best = i_curr;

        i_prev = i_curr;
        i_curr = i_next;
    }

    // Reject split points too close to boundaries.
    if split_point != SPLIT_POINT_INVALID {
        let dist_from_l = if split_point >= index_l {
            split_point - index_l
        } else {
            split_point + points_len - index_l
        };
        let dist_from_r = if index_r >= split_point {
            index_r - split_point
        } else {
            index_r + points_len - split_point
        };
        if dist_from_l < 3 || dist_from_r < 3 {
            split_point = SPLIT_POINT_INVALID;
        }
    }

    split_point
}
