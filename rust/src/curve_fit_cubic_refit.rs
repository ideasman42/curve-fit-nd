///
/// Curve Re-fitting Method
/// =======================
///
/// This is a more processor-intensive method of fitting compared to
/// direct curve fitting, and works as follows:
///
/// - First iteratively remove all points under the error threshold.
///
/// - If corner calculation is enabled:
///   - Find adjacent knots that exceed the angle limit.
///   - Collapse the knots into one, or remove entirely
///     (depending on their error values).
///
/// - Run a re-fit pass, where knots are re-positioned between their adjacent knots
///   when their re-fit position has a lower 'error'.
///   While re-fitting, remove knots that fall below the error threshold.
///
/// This module contains three refinement algorithms:
/// - `refine_remove`: Incremental point removal based on error threshold
/// - `refine_refit`: Re-positioning knots to minimize error
/// - `refine_corner`: Corner detection and collapse
///

use crate::math_vector::{
    dot_vnvn,
    len_squared_vn,
    normalize_vn,
    project_plane_vnvn_normalized,
    sub_vnvn,
};

use super::curve_fit_cubic;

/// Sentinel value representing an invalid index (used for prev/next in non-cyclic curves).
pub const INVALID: usize = ::std::usize::MAX;

/// Type definitions for curve refinement.
pub mod types {
    /// A knot in the curve fitting linked list.
    ///
    /// Knots form a doubly-linked list where each knot represents a control point
    /// in the simplified curve. During refinement, knots can be removed, repositioned,
    /// or marked as corners.
    pub struct Knot {
        /// Index of next knot in the linked list (INVALID if none).
        pub next: usize,
        /// Index of previous knot in the linked list (INVALID if none).
        pub prev: usize,

        /// The index of this knot in the point array.
        ///
        /// Currently the same, access as different for now,
        /// since we may want to support different point/knot indices
        pub index: usize,

        /// If true, this knot cannot be removed (e.g., endpoints of non-cyclic curves).
        pub no_remove: bool,
        /// If true, this knot has been removed from the active list.
        pub is_remove: bool,
        /// If true, this knot is a corner (tangents are not continuous).
        pub is_corner: bool,

        /// Handle lengths [left, right] for the bezier curves meeting at this knot.
        ///
        /// These are signed values: positive extends in the tangent direction,
        /// negative extends opposite. Used to compute control points as:
        /// `control_point = knot_position + tangent * handle_length`
        pub handles: [f64; 2],

        /// Store the error value, to see if we can improve on it
        /// (without having to re-calculate each time)
        ///
        /// This is the error between this knot and the next.
        pub fit_error_sq_next: f64,

        /// Indices into the tangent array [incoming, outgoing].
        ///
        /// Initially these point to contiguous memory (knot_index * 2),
        /// but may be reassigned when corners are created.
        pub tan: [usize; 2],
    }

    /// Immutable point data used during curve refinement.
    ///
    /// Contains references to the original points, cached segment lengths,
    /// and tangent vectors. For cyclic curves, the arrays may be doubled
    /// to allow extracting contiguous slices across the start/end boundary.
    pub struct PointData<'a> {
        /// The input points as a flat array (may be doubled for cyclic curves).
        /// Layout: [p0_x, p0_y, ..., p1_x, p1_y, ..., ...]
        pub points: &'a [f64],
        /// Number of dimensions per point.
        pub dims: usize,
        /// The actual number of unique points.
        pub points_len: usize,

        /// Cached segment lengths between consecutive points.
        /// This array may be doubled as well for cyclic curves.
        pub points_length_cache: &'a [f64],

        /// Tangent vectors at each knot as a flat array (2 per knot: incoming and outgoing).
        /// Layout: [tan0_in_x, tan0_in_y, ..., tan0_out_x, tan0_out_y, ..., ...]
        pub tangents: &'a [f64],
    }

    /// Handle lengths and error for the 2 segments between 3 knots.
    #[derive(Copy, Clone)]
    pub struct KnotAdjacentParams {
        /// Handle lengths [left, right] for the segment before this knot.
        pub handles_prev: [f64; 2],
        /// Handle lengths [left, right] for the segment after this knot.
        pub handles_next: [f64; 2],
        /// Squared error for the segment before this knot.
        pub error_sq_prev: f64,
        /// Squared error for the segment after this knot.
        pub error_sq_next: f64,
    }
}

/// Re-export core types from the types module.
pub use self::types::{
    Knot,
    PointData,
    KnotAdjacentParams,
};

/// Get a point slice from the flat points array.
#[inline]
fn get_point<'a>(pd: &'a PointData, index: usize) -> &'a [f64] {
    let start = index * pd.dims;
    &pd.points[start..start + pd.dims]
}

/// Get a tangent slice from the flat tangents array.
#[inline]
fn get_tangent<'a>(pd: &'a PointData, tan_index: usize) -> &'a [f64] {
    let start = tan_index * pd.dims;
    &pd.tangents[start..start + pd.dims]
}

/// Advance knot index to next knot in array, wrapping at end.
#[inline]
fn knot_step_next_wrap(k_step: &mut usize, knots_end: usize) {
    if *k_step != knots_end {
        *k_step += 1;
    } else {
        // Wrap around.
        *k_step = 0;
    }
}

/// Find the knot furthest from the line between `k_prev` and `k_next`.
///
/// This is used to find a good split point when a curve segment
/// exceeds the error threshold. The point with maximum perpendicular
/// distance from the chord is chosen as the split point.
pub fn knot_find_split_point(
    pd: &PointData,
    knots: &[Knot],
    k_prev: &Knot,
    k_next: &Knot,
) -> usize {
    let mut split_point: usize = INVALID;
    let mut split_point_dist_best: f64 = -::std::f64::MAX;

    let offset = get_point(pd, k_prev.index);

    let mut v_plane = sub_vnvn(get_point(pd, k_prev.index), get_point(pd, k_next.index), pd.dims);
    normalize_vn(&mut v_plane[..pd.dims]);

    let knots_end = knots.len() - 1;
    let mut k_step = k_prev.index;
    loop {
        knot_step_next_wrap(&mut k_step, knots_end);

        if k_step != k_next.index {
            let knot = &knots[k_step];
            let v_offset = sub_vnvn(get_point(pd, knot.index), offset, pd.dims);
            let v_proj = project_plane_vnvn_normalized(&v_offset[..pd.dims], &v_plane[..pd.dims], pd.dims);
            let split_point_dist_test = len_squared_vn(&v_proj[..pd.dims]);
            if split_point_dist_test > split_point_dist_best {
                split_point_dist_best = split_point_dist_test;
                split_point = knot.index;
            }
        } else {
            break;
        }
    }

    split_point
}

/// Find the knot furthest from the line between `k_prev` and `k_next` along a given axis.
///
/// Similar to #knot_find_split_point, but projects points onto the given
/// plane normal instead of perpendicular to the chord. Used for corner
/// detection to find split points that best separate angled segments.
pub fn knot_find_split_point_on_axis(
    pd: &PointData,
    knots: &[Knot],
    k_prev: &Knot,
    k_next: &Knot,
    plane_no: &[f64],
) -> usize {
    let mut split_point: usize = INVALID;
    let mut split_point_dist_best: f64 = -::std::f64::MAX;

    let knots_end = knots.len() - 1;
    let mut k_step = k_prev.index;
    loop {
        knot_step_next_wrap(&mut k_step, knots_end);

        if k_step != k_next.index {
            let knot = &knots[k_step];
            let split_point_dist_test = dot_vnvn(plane_no, get_point(pd, knot.index));
            if split_point_dist_test > split_point_dist_best {
                split_point_dist_best = split_point_dist_test;
                split_point = knot.index;
            }
        } else {
            break;
        }
    }

    split_point
}

/// Find the split point based on sign change of perpendicular distance.
///
/// This finds where the curve crosses the line between the two knots,
/// selecting the crossing point with the largest perpendicular distance.
///
/// For N-dimensional support, we establish a reference perpendicular direction
/// from the first point that deviates from the line, then measure signed distance
/// as the dot product with that reference. This gives consistent "sides" in any dimension.
///
/// Wrapper for #split_point_find_sign_change that uses Knot structures.
pub fn knot_find_split_point_sign_change(
    pd: &PointData,
    _knots: &[Knot],
    k_prev: &Knot,
    k_next: &Knot,
) -> usize {
    curve_fit_cubic::split_point_find_sign_change(
        pd.points,
        pd.points_len,
        k_prev.index,
        k_next.index,
        pd.dims,
    )
}

/// Find the split point with maximum perpendicular distance from the line-segment.
///
/// Wrapper for #split_point_find_max_distance that uses Knot structures.
pub fn knot_find_split_point_max_distance(
    pd: &PointData,
    _knots: &[Knot],
    k_prev: &Knot,
    k_next: &Knot,
) -> usize {
    curve_fit_cubic::split_point_find_max_distance(
        pd.points,
        pd.points_len,
        k_prev.index,
        k_next.index,
        pd.dims,
    )
}

/// Find inflection point where curvature changes sign.
///
/// Wrapper for #split_point_find_inflection that uses Knot structures.
pub fn knot_find_split_point_inflection(
    pd: &PointData,
    _knots: &[Knot],
    k_prev: &Knot,
    k_next: &Knot,
) -> usize {
    curve_fit_cubic::split_point_find_inflection(
        pd.points,
        pd.points_len,
        k_prev.index,
        k_next.index,
        pd.dims,
    )
}

/// Methods for calculating split points during refit.
#[derive(Clone, Copy)]
#[repr(u8)]
enum SplitCalcMethod {
    /// Find the point with maximum error from the fitted curve.
    /// First to try: zero cost (reuses already-calculated refit_index).
    MaxError = 0,
    /// Find the point with maximum perpendicular distance from line-segment.
    /// Good early candidate: always returns a valid result, good general-purpose fallback.
    MaxDistance = 1,
    /// Find inflection point where the curve changes from bending one way to the other.
    /// Useful for S-curves: detects where bending direction changes.
    ///
    /// Try later: may return INVALID.
    Inflection = 2,
    /// Find the point where the curve crosses the line between endpoints (sign change).
    /// Useful for S-curves: detects where the curve crosses the line-segment.
    ///
    /// Try last: may return INVALID.
    SignChange = 3,
}

/// Array of split calculation methods to try, in order.
const SPLIT_CALC_METHODS: &[SplitCalcMethod] = &[
    SplitCalcMethod::MaxError,
    SplitCalcMethod::MaxDistance,
    SplitCalcMethod::Inflection,
    SplitCalcMethod::SignChange,
];

/// Number of split calculation methods.
const SPLIT_CALC_METHODS_NUM: usize = SPLIT_CALC_METHODS.len();


/// Fit a curve segment and return error metrics.
///
/// Returns (error_sq, error_index, [handle_left, handle_right]).
fn knot_remove_error_value(
    dims: usize,
    tan_l: &[f64],
    tan_r: &[f64],
    points_offset: &[f64],
    points_offset_length_cache: &[f64],
) -> (f64, usize, [f64; 2]) {
    let points_len = points_offset.len() / dims;
    let p0 = &points_offset[0..dims];
    let p_last = &points_offset[(points_len - 1) * dims..];

    let ((error_sq, error_index), handle_factor_l, handle_factor_r, _threshold_met) =
        curve_fit_cubic::curve_fit_cubic_to_points_single(
            points_offset, dims, points_offset_length_cache,
            tan_l, tan_r,
            0.0,  // Pass 0.0 to disable early exits (like C does).
            );
    let sub_l = sub_vnvn(&handle_factor_l[..dims], p0, dims);
    let sub_r = sub_vnvn(&handle_factor_r[..dims], p_last, dims);
    (
        error_sq, error_index,
        [dot_vnvn(tan_l, &sub_l[..dims]),
         dot_vnvn(tan_r, &sub_r[..dims])],
    )
}

/// Calculate the number of points from `index_l` to `index_r` inclusive.
///
/// Handles cyclic wrap-around when `index_l > index_r`.
#[inline]
fn knot_span_length(index_l: usize, index_r: usize, points_len: usize) -> usize {
    (if index_l <= index_r {
        index_r - index_l
    } else {
        (index_r + points_len) - index_l
    }) + 1
}

/// Calculate the curve fit error between two knots, including the error index.
///
/// Returns (error_sq, error_index, [handle_left, handle_right]).
/// The error_index is the global index of the point with maximum error.
pub fn knot_calc_curve_error_value_and_index(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64],
    tan_r: &[f64],
) -> (f64, usize, [f64; 2]) {
    let points_offset_len = knot_span_length(knot_l.index, knot_r.index, pd.points_len);

    if points_offset_len != 2 {
        let points_offset_start = knot_l.index * pd.dims;
        let points_offset_end = points_offset_start + points_offset_len * pd.dims;
        let mut result = knot_remove_error_value(
            pd.dims,
            tan_l, tan_r,
            &pd.points[points_offset_start..points_offset_end],
            &pd.points_length_cache[knot_l.index..knot_l.index + points_offset_len],
            );

        // Adjust the offset index to the global index & wrap if needed.
        result.1 += knot_l.index;
        if result.1 >= pd.points_len {
            result.1 -= pd.points_len;
        }

        result
    } else {
        // No points between, use 1/3 handle length with no error as a fallback.
        debug_assert!(points_offset_len == 2);
        let handle_len = pd.points_length_cache[knot_l.index] / 3.0;
        (0.0, knot_l.index, [handle_len, handle_len])
    }
}

/// Calculate the curve fit error between two knots (without error index).
///
/// Returns (error_sq, [handle_left, handle_right]).
pub fn knot_calc_curve_error_value(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64],
    tan_r: &[f64],
) -> (f64, [f64; 2]) {
    let points_offset_len = knot_span_length(knot_l.index, knot_r.index, pd.points_len);

    if points_offset_len != 2 {
        let points_offset_start = knot_l.index * pd.dims;
        let points_offset_end = points_offset_start + points_offset_len * pd.dims;
        let result = knot_remove_error_value(
            pd.dims,
            tan_l, tan_r,
            &pd.points[points_offset_start..points_offset_end],
            &pd.points_length_cache[knot_l.index..knot_l.index + points_offset_len],
            );
        (result.0, result.2)
    } else {
        // No points between, use 1/3 handle length with no error as a fallback.
        debug_assert!(points_offset_len == 2);
        let handle_len = pd.points_length_cache[knot_l.index] / 3.0;
        (0.0, [handle_len, handle_len])
    }
}

macro_rules! unlikely { ($body:expr) => { $body } }

/// First refinement pass: iteratively remove knots below the error threshold.
pub mod refine_remove {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
        get_tangent,
        Knot,
        PointData,
    };
    use crate::min_heap;

    /// State stored in the heap for potential knot removal.
    ///
    /// Stores the handle lengths that would result from removing this knot.
    #[derive(Copy, Clone)]
    struct KnotRemoveState {
        /// Index of the knot being considered for removal.
        index: usize,
        /// Handle lengths if this knot is removed.
        handles: [f64; 2],
    }

    /// (Re)calculate the error for removing a knot and update the heap.
    ///
    /// Tests the error that would result from removing k_curr and fitting
    /// a single curve from k_prev to k_next. Updates the heap if the error
    /// is below the threshold.
    fn knot_remove_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<(f64, usize), KnotRemoveState>,
        knots: &[Knot],
        knots_handle: &mut [min_heap::NodeHandle],
        k_curr: &Knot,
        error_max_sq: f64,
    ) {
        debug_assert!(k_curr.no_remove == false);

        let (fit_error_max_sq, handles) = {
            let k_prev = &knots[k_curr.prev];
            let k_next = &knots[k_curr.next];

            let result = knot_calc_curve_error_value(
                pd, k_prev, k_next,
                get_tangent(pd, k_prev.tan[1]),
                get_tangent(pd, k_next.tan[0]));
            result
        };

        let k_curr_heap_node = &mut knots_handle[k_curr.index];
        if fit_error_max_sq < error_max_sq {
            heap.insert_or_update(
                k_curr_heap_node,
                (fit_error_max_sq, k_curr.index),
                KnotRemoveState {
                    index: k_curr.index,
                    handles: handles,
                },
            );
        } else {
            if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                heap.remove(*k_curr_heap_node);
                *k_curr_heap_node = min_heap::NodeHandle::INVALID;
            }
        }
    }

    /// Iteratively remove all points under the error threshold.
    ///
    /// Uses a min-heap to efficiently find and remove the knot with the
    /// smallest error value. After each removal, the error values of
    /// adjacent knots are recalculated and the heap is updated.
    ///
    /// This is the first pass of the curve refinement algorithm.
    pub fn curve_incremental_simplify(
        pd: &PointData,
        knots: &mut [Knot],
        knots_handle: &mut [min_heap::NodeHandle],
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
    ) {
        // Use (error, index) tuple for deterministic ordering when errors are equal.
        // This matches C's implicit ordering from insertion order in the binary heap.
        let mut heap = min_heap::MinHeap::<(f64, usize), KnotRemoveState>::with_capacity(knots.len());

        for k_index in 0..knots.len() {
            let k_curr = &knots[k_index];
            if (k_curr.no_remove == false) &&
               (k_curr.is_remove == false) &&
               (k_curr.is_corner == false)
            {
                knot_remove_error_recalculate(
                    pd, &mut heap, knots, knots_handle, k_curr, error_max_sq);
            }
        }

        while let Some(((error_sq, _), r)) = heap.pop_min_with_value() {
            knots_handle[r.index] = min_heap::NodeHandle::INVALID;

            let k_next_index;
            let k_prev_index;
            {
                // let r: &mut remove_states[r_index];
                let k_curr: &mut Knot = &mut knots[r.index];

                if unlikely!(*knots_len_remaining <= 2) {
                    continue;
                }

                k_next_index = k_curr.next;
                k_prev_index = k_curr.prev;

                k_curr.is_remove = true;

                if cfg!(debug_assertions) {
                    k_curr.next = INVALID;
                    k_curr.prev = INVALID;
                }
            }
            knots[k_prev_index].handles[1] = r.handles[0];
            knots[k_next_index].handles[0] = r.handles[1];

            debug_assert!(error_sq <= error_max_sq);

            knots[k_prev_index].fit_error_sq_next = error_sq;
            // Remove ourselves.
            knots[k_next_index].prev = k_prev_index;
            knots[k_prev_index].next = k_next_index;


            for k_iter_index in &[k_prev_index, k_next_index] {
                let k_iter = &knots[*k_iter_index];
                if (k_iter.no_remove == false) &&
                   (k_iter.is_corner == false) &&
                   (k_iter.prev != INVALID) &&
                   (k_iter.next != INVALID)
                {
                    knot_remove_error_recalculate(
                        pd, &mut heap, knots, knots_handle, k_iter, error_max_sq);
                }
            }

            *knots_len_remaining -= 1;
        }
        drop(heap);
    }
}
// end refine_remove


/// Second refinement pass: reposition knots to minimize error.
///
/// Tests moving each knot to positions between its neighbors to find
/// optimal placement. Knots may also be removed if they fall below threshold.
pub mod refine_refit {

    use super::{
        INVALID,
        knot_calc_curve_error_value,
        knot_calc_curve_error_value_and_index,
        knot_find_split_point,
        get_tangent,
        Knot,
        KnotAdjacentParams,
        PointData,
    };
    use crate::min_heap;

    /// Result from refining a refit index in one direction.
    struct RefineResult {
        /// The best refit index found during the search.
        index_refit: usize,
        /// Handle lengths and errors for the refit position.
        params: KnotAdjacentParams,
        /// True if an improvement was found over the initial position.
        is_refined: bool,
    }

    /// Refine the refit index by searching neighbors for lower error.
    /// Stops when no further improvement is found.
    /// `dir`: -1 to search toward `k_prev`, 1 to search toward `k_next`.
    fn knot_refit_index_refine(
        pd: &PointData,
        knots: &[Knot],
        k_prev: &Knot,
        k_next: &Knot,
        index_refit: usize,
        dir: isize,
        mut cost_sq_max: f64,
    ) -> RefineResult {
        // Stop before reaching the adjacent knot.
        let index_end = if dir == -1 { k_prev.index } else { k_next.index };
        let points_len = pd.points_len;
        let mut i = index_refit;
        let mut result = RefineResult {
            index_refit,
            params: KnotAdjacentParams {
                handles_prev: [0.0; 2],
                handles_next: [0.0; 2],
                error_sq_prev: 0.0,
                error_sq_next: 0.0,
            },
            is_refined: false,
        };

        // Step through indices in direction `dir`, with wraparound.
        loop {
            i = ((i as isize + dir + points_len as isize) as usize) % points_len;
            if i == index_end {
                break;
            }

            let k_test = &knots[i];
            let (error_sq_prev, handles_prev_test) = knot_calc_curve_error_value(
                pd, k_prev, k_test,
                get_tangent(pd, k_prev.tan[1]),
                get_tangent(pd, k_test.tan[0]),
            );
            if error_sq_prev >= cost_sq_max {
                break;
            }

            let (error_sq_next, handles_next_test) = knot_calc_curve_error_value(
                pd, k_test, k_next,
                get_tangent(pd, k_test.tan[1]),
                get_tangent(pd, k_next.tan[0]),
            );
            if error_sq_next >= cost_sq_max {
                break;
            }

            // Raise the bar: subsequent iterations must beat this.
            cost_sq_max = error_sq_prev.max(error_sq_next);
            result.index_refit = i;
            result.params.handles_prev = handles_prev_test;
            result.params.handles_next = handles_next_test;
            result.params.error_sq_prev = error_sq_prev;
            result.params.error_sq_next = error_sq_next;
            result.is_refined = true;
        }
        result
    }

    /// State stored in the heap for potential knot repositioning.
    #[derive(Copy, Clone)]
    struct KnotRefitState {
        /// Index of the knot being considered for refit.
        index: usize,
        /// Target position index for refitting. When INVALID, remove this knot instead.
        index_refit: usize,
        /// Handle lengths and errors for the refit position.
        fit_params: KnotAdjacentParams,
    }

    /// (Re)calculate the error and optimal refit position for a knot.
    fn knot_refit_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<(f64, usize), KnotRefitState>,
        knots: &[Knot],
        knots_handle: &mut [min_heap::NodeHandle],
        k_curr: &Knot,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
        use_refit_remove: bool,
    ) {
        debug_assert!(k_curr.no_remove == false);

        let k_curr_heap_node = &mut knots_handle[k_curr.index];

        let k_prev = &knots[k_curr.prev];
        let k_next = &knots[k_curr.next];

        let mut k_refit_index;

        if use_refit_remove {
            // Support re-fitting to remove points.
            let (fit_error_max_sq, fit_error_index, handles) =
                knot_calc_curve_error_value_and_index(
                    pd, k_prev, k_next,
                    get_tangent(pd, k_prev.tan[1]),
                    get_tangent(pd, k_next.tan[0]),
                    );

            if fit_error_max_sq < error_max_sq {
                // Always perform removal before refitting, (make a negative number)
                heap.insert_or_update(
                    k_curr_heap_node,
                    // Weight for the greatest improvement, with index as tiebreaker.
                    (fit_error_max_sq - error_max_sq, k_curr.index),
                    KnotRefitState {
                        index: k_curr.index,
                        // INVALID == remove
                        index_refit: INVALID,
                        fit_params: KnotAdjacentParams {
                            handles_prev: [handles[0], 0.0],  // [1] unused
                            handles_next: [0.0, handles[1]],  // [0] unused
                            error_sq_prev: fit_error_max_sq,
                            error_sq_next: fit_error_max_sq,
                        },
                    }
                );
                return;
            }

            // Use the largest point of difference when removing
            // as the target to refit to.
            k_refit_index = fit_error_index;
        } else {
            k_refit_index = knot_find_split_point(pd, knots, k_prev, k_next);
        }

        let cost_sq_src_max = k_prev.fit_error_sq_next.max(k_curr.fit_error_sq_next);
        debug_assert!(cost_sq_src_max <= error_max_sq);

        /// Calculate curve errors for both segments around a refit position.
        /// Returns None if either segment exceeds the error threshold.
        fn knot_calc_curve_error_value_pair_above_error_or_none(
            pd: &PointData, k_prev: &Knot, k_refit: &Knot, k_next: &Knot, error_max_sq: f64,
        ) -> Option<([f64; 2], f64, [f64; 2], f64)> {
            let (fit_error_prev, handles_prev) =
                knot_calc_curve_error_value(
                    pd, k_prev, k_refit,
                    get_tangent(pd, k_prev.tan[1]),
                    get_tangent(pd, k_refit.tan[0]),
                );

            if fit_error_prev < error_max_sq {
                let (fit_error_next, handles_next) =
                    knot_calc_curve_error_value(
                        pd, k_refit, k_next,
                        get_tangent(pd, k_refit.tan[1]),
                        get_tangent(pd, k_next.tan[0]),
                    );
                if fit_error_next < error_max_sq {
                    return Some((
                        handles_prev, fit_error_prev,
                        handles_next, fit_error_next,
                    ));
                }
            }
            None
        }

        // cache result of 'knot_calc_curve_error_value_pair_above_error_or_none'
        let mut refit_result_or_none: Option<([f64; 2], f64, [f64; 2], f64)> = None;

        if use_optimize_exhaustive {

            // loop over inner knots
            let mut k_test_index = k_prev.index + 1;

            // start with current state
            let mut cost_sq_best = cost_sq_src_max;

            loop {
                if k_test_index == knots.len() {
                    k_test_index = 0;
                }
                if k_test_index == k_next.index {
                    break;
                }

                if k_test_index != k_curr.index {
                    if let Some(fit_result_test) =
                        knot_calc_curve_error_value_pair_above_error_or_none(
                            pd, k_prev, &knots[k_test_index], k_next, cost_sq_best)
                    {
                        let cost_sq_test_prev = fit_result_test.1;
                        let cost_sq_test_next = fit_result_test.3;
                        cost_sq_best = cost_sq_test_prev.max(cost_sq_test_next);
                        k_refit_index = k_test_index;

                        // Result for re-use if this is the best fit.
                        refit_result_or_none = Some(fit_result_test);
                    }
                }
                k_test_index += 1;
            }
        } else {
            // Try multiple split calculation methods and pick the best one.
            let mut best_cost_sq_max = f64::MAX;

            // Track all indices tried to avoid redundant error calculations.
            let mut tried_indices: [usize; super::SPLIT_CALC_METHODS_NUM] = [INVALID; super::SPLIT_CALC_METHODS_NUM];
            let mut tried_indices_num: usize = 0;

            for method in super::SPLIT_CALC_METHODS {
                let test_refit_index = match method {
                    super::SplitCalcMethod::MaxError => k_refit_index,  // Already calculated above
                    super::SplitCalcMethod::SignChange => super::knot_find_split_point_sign_change(pd, knots, k_prev, k_next),
                    super::SplitCalcMethod::MaxDistance => super::knot_find_split_point_max_distance(pd, knots, k_prev, k_next),
                    super::SplitCalcMethod::Inflection => super::knot_find_split_point_inflection(pd, knots, k_prev, k_next),
                };

                if test_refit_index == INVALID || test_refit_index == k_curr.index {
                    continue;
                }

                // Skip if this index was already evaluated by a previous method.
                let already_tried = tried_indices[..tried_indices_num].contains(&test_refit_index);
                if already_tried {
                    continue;
                }
                tried_indices[tried_indices_num] = test_refit_index;
                tried_indices_num += 1;

                if let Some(fit_result_test) =
                    knot_calc_curve_error_value_pair_above_error_or_none(
                        pd, k_prev, &knots[test_refit_index], k_next, cost_sq_src_max)
                {
                    let test_cost_sq_max = fit_result_test.1.max(fit_result_test.3);
                    if test_cost_sq_max < best_cost_sq_max {
                        best_cost_sq_max = test_cost_sq_max;
                        k_refit_index = test_refit_index;
                        refit_result_or_none = Some(fit_result_test);

                        // Perfect fit, no point trying other methods.
                        if best_cost_sq_max == 0.0 {
                            break;
                        }
                    }
                }
            }

            // Local refinement: search neighbors for a better refit index.
            if let Some((_handles_prev, fit_error_dst_prev, _handles_next, fit_error_dst_next)) =
                refit_result_or_none
            {
                let cost_sq_dst_max_init = fit_error_dst_prev.max(fit_error_dst_next);

                if cost_sq_dst_max_init > 0.0 {
                    // Search toward k_prev (dir=-1) and toward k_next (dir=1).
                    let scan_prev = knot_refit_index_refine(
                        pd, knots, k_prev, k_next, k_refit_index, -1, cost_sq_dst_max_init);
                    let scan_next = knot_refit_index_refine(
                        pd, knots, k_prev, k_next, k_refit_index, 1, cost_sq_dst_max_init);

                    // Pick the best result from both directions.
                    if scan_prev.is_refined || scan_next.is_refined {
                        let scan = if scan_prev.is_refined && scan_next.is_refined {
                            let cost_sq_max_prev = scan_prev.params.error_sq_prev.max(scan_prev.params.error_sq_next);
                            let cost_sq_max_next = scan_next.params.error_sq_prev.max(scan_next.params.error_sq_next);

                            if cost_sq_max_prev < cost_sq_max_next {
                                &scan_prev
                            } else if cost_sq_max_next < cost_sq_max_prev {
                                &scan_next
                            } else {
                                let cost_sq_min_prev = scan_prev.params.error_sq_prev.min(scan_prev.params.error_sq_next);
                                let cost_sq_min_next = scan_next.params.error_sq_prev.min(scan_next.params.error_sq_next);
                                if cost_sq_min_prev <= cost_sq_min_next { &scan_prev } else { &scan_next }
                            }
                        } else if scan_prev.is_refined {
                            &scan_prev
                        } else {
                            &scan_next
                        };

                        // Use results from the winning direction.
                        k_refit_index = scan.index_refit;
                        refit_result_or_none = Some((
                            scan.params.handles_prev, scan.params.error_sq_prev,
                            scan.params.handles_next, scan.params.error_sq_next,
                        ));
                    }
                }
            }
        }
        // end exhaustive test

        if let Some((
            handles_prev, fit_error_dst_prev,
            handles_next, fit_error_dst_next,
        )) = refit_result_or_none {
            let fit_error_dst_max_sq =
                fit_error_dst_prev.max(fit_error_dst_next);
            debug_assert!(fit_error_dst_max_sq < cost_sq_src_max);
            heap.insert_or_update(
                k_curr_heap_node,
                // Weight for the greatest improvement, with index as tiebreaker.
                (cost_sq_src_max - fit_error_dst_max_sq, k_curr.index),
                KnotRefitState {
                    index: k_curr.index,
                    index_refit: k_refit_index,
                    fit_params: KnotAdjacentParams {
                        handles_prev: handles_prev,
                        handles_next: handles_next,
                        error_sq_prev: fit_error_dst_prev,
                        error_sq_next: fit_error_dst_next,
                    },
                }
            );
            return;
        }

        if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
            heap.remove(*k_curr_heap_node);
            *k_curr_heap_node = min_heap::NodeHandle::INVALID;
        }
    }

    /// Re-adjust the curves by re-fitting points.
    pub fn curve_incremental_simplify_refit(
        pd: &PointData,
        knots: &mut [Knot],
        knots_handle: &mut [min_heap::NodeHandle],
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
        use_refit_remove: bool,
    ) {
        // Use (error, index) tuple for deterministic ordering when errors are equal.
        let mut heap =
            min_heap::MinHeap::<(f64, usize), KnotRefitState>::with_capacity(*knots_len_remaining);

        let mut _added_count = 0;
        for k_index in 0..knots.len() {
            let k_curr = &knots[k_index];
            if (k_curr.no_remove == false) &&
               (k_curr.is_remove == false) &&
               (k_curr.is_corner == false)
            {
                let old_handle = knots_handle[k_index];
                knot_refit_error_recalculate(
                    pd, &mut heap, knots, knots_handle, k_curr,
                    error_max_sq, use_optimize_exhaustive, use_refit_remove);
                if knots_handle[k_index] != old_handle {
                    _added_count += 1;
                }
            }
        }

        while let Some(r) = heap.pop_min() {
            knots_handle[r.index] = min_heap::NodeHandle::INVALID;

            let k_prev_index;
            let k_next_index;
            {
                {
                    let k_old = &knots[r.index];
                    k_prev_index = k_old.prev;
                    k_next_index = k_old.next;
                }

                knots[k_prev_index].handles[1] = r.fit_params.handles_prev[0];
                knots[k_next_index].handles[0] = r.fit_params.handles_next[1];

                // Update error values for changed segments.
                //
                // Before:
                // - `k_prev - (error_sq_prev) -> k_refit - (error_sq_next) -> k_next`.
                // After:
                // - `k_prev->fit_error_sq_next := error_sq_prev`.
                // - `k_refit->fit_error_sq_next := error_sq_next`.
                // - `k_next->fit_error_sq_next`: unchanged (segment beyond k_next unaffected).
                knots[k_prev_index].fit_error_sq_next = r.fit_params.error_sq_prev;

                if r.index_refit != INVALID {
                    let k_refit = &mut knots[r.index_refit];
                    k_refit.handles[0] = r.fit_params.handles_prev[1];
                    k_refit.handles[1] = r.fit_params.handles_next[0];
                    k_refit.fit_error_sq_next = r.fit_params.error_sq_next;
                }
            }
            // finished with 'r'

            // Skip if curve is too small to simplify further.
            if unlikely!(*knots_len_remaining <= 2) {
                continue;
            }

            {
                let k_old = &mut knots[r.index];
                k_old.next = INVALID;
                k_old.prev = INVALID;
                k_old.is_remove = true;
            }

            if r.index_refit == INVALID {
                knots[k_next_index].prev = k_prev_index;
                knots[k_prev_index].next = k_next_index;

                *knots_len_remaining -= 1;
            } else {
                // Remove ourselves.
                knots[k_next_index].prev = r.index_refit;
                knots[k_prev_index].next = r.index_refit;

                let k_refit = &mut knots[r.index_refit];
                k_refit.prev = k_prev_index;
                k_refit.next = k_next_index;

                k_refit.is_remove = false;
            }

            for k_iter_index in &[k_prev_index, k_next_index] {
                let k_iter = &knots[*k_iter_index];
                if (k_iter.no_remove == false) &&
                   (k_iter.is_corner == false) &&
                   (k_iter.prev != INVALID) &&
                   (k_iter.next != INVALID)
                {
                    knot_refit_error_recalculate(
                        pd, &mut heap, knots, knots_handle, k_iter,
                        error_max_sq, use_optimize_exhaustive, use_refit_remove);
                }
            }
        }

        drop(heap);
    }
}
// end refine_refit

/// Corner detection pass: identify and collapse sharp angle transitions.
pub mod refine_corner {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
        knot_find_split_point_on_axis,
        get_point,
        get_tangent,
        Knot,
        KnotAdjacentParams,
        PointData,
    };
    use crate::math_vector::{
        dot_vnvn,
        len_squared_vnvn,
        project_vnvn_normalized,
        sub_vnvn,
    };
    use crate::min_heap;

    /// Result of collapsing a corner.
    #[derive(Copy, Clone)]
    struct KnotCornerState {
        /// Index of the knot being considered for corner collapse.
        index: usize,
        /// Indices of adjacent knots [k_prev, k_next] whose tangents will be
        /// collapsed into this corner.
        index_pair: [usize; 2],
        /// Handle lengths and errors for the collapsed corner.
        fit_params: KnotAdjacentParams,
    }

    /// (Re)calculate the error incurred from turning this into a corner.
    fn knot_corner_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<(f64, usize), KnotCornerState>,
        knots_handle: &mut [min_heap::NodeHandle],
        k_split: &Knot,
        k_prev: &Knot,
        k_next: &Knot,
        error_max_sq: f64,
    ) {
        debug_assert!(
            (k_prev.no_remove == false) &&
            (k_next.no_remove == false)
        );

        let k_split_heap_node = &mut knots_handle[k_split.index];

        // Test skipping 'k_prev' by using points (k_prev.prev to k_split).
        {
            let (fit_error_dst_prev, handles_prev) =
                knot_calc_curve_error_value(
                    pd, k_prev, k_split,
                    get_tangent(pd, k_prev.tan[1]),
                    get_tangent(pd, k_prev.tan[1]),
                    );
            if fit_error_dst_prev < error_max_sq {
                let (fit_error_dst_next, handles_next) =
                    knot_calc_curve_error_value(
                        pd, k_split, k_next,
                        get_tangent(pd, k_next.tan[0]),
                        get_tangent(pd, k_next.tan[0]),
                        );
                if fit_error_dst_next < error_max_sq {
                    // _must_ be assigned to k_split, later
                    heap.insert_or_update(
                        k_split_heap_node,
                        // Weight for the greatest improvement, with index as tiebreaker.
                        (fit_error_dst_prev.max(fit_error_dst_next), k_split.index),
                        KnotCornerState {
                            index: k_split.index,
                            index_pair: [k_prev.index, k_next.index],
                            fit_params: KnotAdjacentParams {
                                handles_prev: handles_prev,
                                handles_next: handles_next,
                                error_sq_prev: fit_error_dst_prev,
                                error_sq_next: fit_error_dst_next,
                            },
                        }
                    );

                    return;
                }
            }
        }

        if *k_split_heap_node != min_heap::NodeHandle::INVALID {
            heap.remove(*k_split_heap_node);
            *k_split_heap_node = min_heap::NodeHandle::INVALID;
        }
    }

    /// Attempt to collapse close knots into corners.
    pub fn curve_incremental_simplify_corners(
        pd: &PointData,
        knots: &mut [Knot],
        knots_handle: &mut [min_heap::NodeHandle],
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
        error_sq_collapse_max: f64,
        corner_angle: f64,
    ) {
        // don't pre-allocate, since its likely there are no corners
        // Use (error, index) tuple for deterministic ordering when errors are equal.
        let mut heap = min_heap::MinHeap::<(f64, usize), KnotCornerState>::with_capacity(0);

        let corner_angle_cos = corner_angle.cos();

        for k_prev_index in 0..knots.len() {
            if let Some((k_prev, k_next)) = {
                let k_prev: &Knot = &knots[k_prev_index];

                if (k_prev.is_remove == false) &&
                   (k_prev.no_remove == false) &&
                   (k_prev.next != INVALID) &&
                   (knots[k_prev.next].no_remove == false)
                {
                    Some((k_prev, &knots[k_prev.next]))
                } else {
                    None
                }
            }
            {
                // Angle outside threshold
                if dot_vnvn(
                    get_tangent(pd, k_prev.tan[0]),
                    get_tangent(pd, k_next.tan[1])) < corner_angle_cos
                {
                    // Measure distance projected onto a plane,
                    //since the points may be offset along their own tangents.
                    let plane_no = sub_vnvn(
                        get_tangent(pd, k_next.tan[0]),
                        get_tangent(pd, k_prev.tan[1]),
                        pd.dims,
                        );

                    // Compare 2x so as to allow both to be changed
                    // by maximum of `error_sq_collapse_max`.
                    let k_split_index = knot_find_split_point_on_axis(
                        pd,
                        knots,
                        k_prev,
                        k_next,
                        &plane_no[..pd.dims],
                        );

                    if k_split_index != INVALID {
                        let co_prev  = get_point(pd, k_prev.index);
                        let co_next  = get_point(pd, k_next.index);
                        let co_split = get_point(pd, k_split_index);

                        let k_proj_ref = project_vnvn_normalized(
                            co_prev, get_tangent(pd, k_prev.tan[1]), pd.dims);
                        let k_proj_split = project_vnvn_normalized(
                            co_split, get_tangent(pd, k_prev.tan[1]), pd.dims);

                        if len_squared_vnvn(
                            &k_proj_ref[..pd.dims], &k_proj_split[..pd.dims]) < error_sq_collapse_max
                        {
                            let k_proj_ref = project_vnvn_normalized(
                                co_next, get_tangent(pd, k_next.tan[0]), pd.dims);
                            let k_proj_split = project_vnvn_normalized(
                                co_split, get_tangent(pd, k_next.tan[0]), pd.dims);

                            if len_squared_vnvn(
                                &k_proj_ref[..pd.dims], &k_proj_split[..pd.dims]) < error_sq_collapse_max
                            {
                                knot_corner_error_recalculate(
                                    pd,
                                    &mut heap,
                                    knots_handle,
                                    &knots[k_split_index],
                                    k_prev,
                                    k_next,
                                    error_max_sq,
                                    );
                            }
                        }
                    }
                }
            }
        }

        while let Some(c) = heap.pop_min() {
            knots_handle[c.index] = min_heap::NodeHandle::INVALID;

            let k_split_index = c.index;
            let k_prev_index = c.index_pair[0];
            let k_next_index = c.index_pair[1];

            let tan_prev;
            let tan_next;


            {
                let k_prev  = &mut knots[k_prev_index];
                k_prev.next = k_split_index;
                k_prev.handles[1]  = c.fit_params.handles_prev[0];
                tan_prev = k_prev.tan[1];

                debug_assert!(c.fit_params.error_sq_prev <= error_max_sq);
                k_prev.fit_error_sq_next = c.fit_params.error_sq_prev;
            }

            {
                let k_next  = &mut knots[k_next_index];
                k_next.prev = k_split_index;
                tan_next = k_next.tan[0];

                k_next.handles[0] = c.fit_params.handles_next[1];
            }

            // Remove while collapsing
            {
                let k_split = &mut knots[k_split_index];

                // Insert
                k_split.is_remove = false;
                k_split.is_corner = true;

                k_split.prev = k_prev_index;
                k_split.next = k_next_index;

                // Update tangents
                k_split.tan[0] = tan_prev;
                k_split.tan[1] = tan_next;

                // Own handles
                k_split.handles[0] = c.fit_params.handles_prev[1];
                k_split.handles[1] = c.fit_params.handles_next[0];

                debug_assert!(c.fit_params.error_sq_next <= error_max_sq);
                k_split.fit_error_sq_next = c.fit_params.error_sq_next;
            }

            *knots_len_remaining += 1;
        }


        drop(heap);
    }
}

// end refine_corner
