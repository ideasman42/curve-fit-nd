
/* TODO, collapse at end-points for non-cyclic curves! - Dont just collapse curves */

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <assert.h>

#include <string.h>
#include <stdlib.h>


#include <stdio.h>

#include "curve_fit_inline.h"
#include "../curve_fit_nd.h"

#include "heap_impl.h"

#ifdef _MSC_VER
#  define alloca(size) _alloca(size)
#endif

#if !defined(_MSC_VER)
#  define USE_VLA
#endif

#ifdef USE_VLA
#  ifdef __GNUC__
#    pragma GCC diagnostic ignored "-Wvla"
#  endif
#else
#  ifdef __GNUC__
#    pragma GCC diagnostic error "-Wvla"
#  endif
#endif







#define USE_CORNER_DETECT

#define SQUARE(a)  ({ \
	typeof(a) a_ = (a); \
	((a_) * (a_)); })

typedef unsigned int uint;


typedef struct Knot {
	struct Knot *next, *prev;

	HeapNode *heap_node;

	uint point_index;  /* index in point array */
	uint knot_index;   /* index in knot array*/
	double handles[2];

	uint can_remove : 1;
	uint is_removed : 1;
	uint is_corner  : 1;

	/* initially point to contiguous memory, however we may re-assign */
	double *tan[2];
} Knot;

struct KnotRemoveState {
	uint knot_index;
	/* handles for prev/next knots */
	double handles[2];
};

static double knot_remove_error_value(
        const double *tan_l, const double *tan_r,
        const double *points_offset, const uint points_offset_len,
        const uint dims,
        /* avoid having to re-calculate again */
        double r_handle_factors[2])
{
	double error_sq = FLT_MAX;
	double handle_factors[2][dims];

	curve_fit_cubic_to_points_single_db(
	        points_offset, points_offset_len, dims, 0.0f,
	        tan_l, tan_r,
	        handle_factors[0], handle_factors[1],
	        &error_sq);

	isub_vnvn(handle_factors[0], points_offset, dims);
	r_handle_factors[0] = dot_vnvn(tan_l, handle_factors[0], dims);

	isub_vnvn(handle_factors[1], &points_offset[(points_offset_len - 1) * dims], dims);
	r_handle_factors[1] = dot_vnvn(tan_r, handle_factors[1], dims);

	return error_sq;
}

static double knot_calc_curve_error_value(
        const double *points, const uint points_len,
        const struct Knot *knot_l, const struct Knot *knot_r,
        const double *tan_l, const double *tan_r,
        const uint dims,
        double r_handle_factors[2])
{
	const double *points_offset;
	uint points_offset_len;

	if (knot_l->point_index < knot_r->point_index) {
		points_offset = &points[knot_l->point_index * dims];
		points_offset_len = (knot_r->point_index - knot_l->point_index) + 1;
	}
	else {
		points_offset = &points[knot_l->point_index * dims];
		points_offset_len = ((knot_r->point_index + points_len) - knot_l->point_index) + 1;
	}

	return knot_remove_error_value(
	        tan_l, tan_r,
	        points_offset, points_offset_len,
	        dims,
	        r_handle_factors);
}

static void knot_remove_error_recalculate(
        Heap *heap, const double *points, const uint points_len,
        struct Knot *k, const double error_sq_max,
        const uint dims)
{
	assert(k->can_remove);
	double handles[2];

	const double cost_sq = knot_calc_curve_error_value(
	        points, points_len,
	        k->prev, k->next,
	        k->prev->tan[1], k->next->tan[0],
	        dims,
	        handles);

	if (cost_sq < error_sq_max) {
		struct KnotRemoveState *r;
		if (k->heap_node) {
			r = HEAP_node_ptr(k->heap_node);
			HEAP_remove(heap, k->heap_node);
		}
		else {
			r = malloc(sizeof(*r));
			r->knot_index = k->knot_index;
		}

		r->handles[0] = handles[0];
		r->handles[1] = handles[1];

		k->heap_node = HEAP_insert(heap, cost_sq, r);
	}
	else {
		if (k->heap_node) {
			struct KnotRemoveState *r;
			r = HEAP_node_ptr(k->heap_node);
			HEAP_remove(heap, k->heap_node);

			free(r);

			k->heap_node = NULL;
		}
	}
}

/**
 * Return length after being reduced.
 */
static uint curve_incremental_simplify(
        const double *points, const uint points_len,
        struct Knot *knots, const uint knots_len, uint knots_len_remaining,
        double error_sq_max, const uint dims)
{
	Heap *heap = HEAP_new(knots_len);
	for (uint i = 0; i < knots_len; i++) {
		struct Knot *k = &knots[i];
		if (k->can_remove && (k->is_removed == false) && (k->is_corner == false)) {
			knot_remove_error_recalculate(heap, points, points_len, k, error_sq_max, dims);
		}
	}

	while (HEAP_is_empty(heap) == false) {
		struct Knot *k;

		{
			struct KnotRemoveState *r = HEAP_popmin(heap);
			k = &knots[r->knot_index];
			k->heap_node = NULL;
			k->prev->handles[1] = r->handles[0];
			k->next->handles[0] = r->handles[1];
			free(r);
		}

		struct Knot *k_prev = k->prev;
		struct Knot *k_next = k->next;

		/* remove ourselves */
		k_next->prev = k_prev;
		k_prev->next = k_next;

		k->next = NULL;
		k->prev = NULL;
		k->is_removed = true;

		if (k_prev->can_remove) {
			knot_remove_error_recalculate(heap, points, points_len, k_prev, error_sq_max, dims);
		}

		if (k_next->can_remove) {
			knot_remove_error_recalculate(heap, points, points_len, k_next, error_sq_max, dims);
		}

		knots_len_remaining -= 1;
	}

	HEAP_free(heap, free);

	return knots_len_remaining;
}

#ifdef USE_CORNER_DETECT

/**
 * Result of collapsing a corner.
 */
struct KnotCornerState {
	uint knot_index;
	/* merge adjacent handles into this one (may be shared with the 'knot_index') */
	uint knot_index_adjacent[2];

	/* handles for prev/next knots */
	double handles_prev[2], handles_next[2];
};

/**
 * Find the knot furthest from the line between \a knot_l & \a knot_r.
 * This is to be used as a split point.
 */
static unsigned int knot_find_split_point(
        const double *points, const uint points_len,
        const struct Knot *knot_l, const struct Knot *knot_r,
        const uint knots_len,
        const double *plane_no,
        const uint dims)
{
	(void)points_len;

	unsigned int split_point = (unsigned int)-1;
	double split_point_dist_best = -DBL_MAX;

	const uint knots_end = knots_len - 1;
	const struct Knot *k_step = knot_l;
	do {
		if (k_step->knot_index != knots_end) {
			k_step += 1;
		}
		else {
			/* wrap around */
			k_step = k_step - knots_end;
		}

		if (k_step != knot_r) {
			double split_point_dist_test = dot_vnvn(plane_no, &points[k_step->knot_index * dims], dims);
			if (split_point_dist_test > split_point_dist_best) {
				split_point_dist_best = split_point_dist_test;
				split_point = k_step->knot_index;
			}
		}
		else {
			break;
		}

	} while (true);

	return split_point;
}

/**
 * (Re)calculate the error incurred from turning this into a corner.
 */
static void knot_corner_error_recalculate(
        Heap *heap, const double *points, const uint points_len,
        struct Knot *k_split,
        struct Knot *k_prev, struct Knot *k_next,
        const double error_sq_max,
        const uint dims)
{
	assert(k_prev->can_remove && k_next->can_remove);

	double handles_prev[2], handles_next[2];

	/* test skipping 'k_prev' by using points (k_prev->prev to k_split) */
	double cost_prev_sq, cost_next_sq;

	if (((cost_prev_sq = knot_calc_curve_error_value(
	          points, points_len,
	          k_prev->prev, k_split,
	          k_prev->prev->tan[1], k_prev->tan[0],
	          dims,
	          handles_prev)) < error_sq_max) &&
	    ((cost_next_sq = knot_calc_curve_error_value(
	          points, points_len,
	          k_split, k_next->next,
	          k_next->tan[1], k_next->next->tan[0],
	          dims,
	          handles_next)) < error_sq_max))
	{
		struct KnotCornerState *c;
		if (k_split->heap_node) {
			c = HEAP_node_ptr(k_split->heap_node);
			HEAP_remove(heap, k_split->heap_node);
		}
		else {
			c = malloc(sizeof(*c));
			c->knot_index = k_split->knot_index;
		}

		c->knot_index_adjacent[0] = k_prev->knot_index;
		c->knot_index_adjacent[1] = k_next->knot_index;

		/* need to store handle lengths for both sides */
		c->handles_prev[0] = handles_prev[0];
		c->handles_prev[1] = handles_prev[1];

		c->handles_next[0] = handles_next[0];
		c->handles_next[1] = handles_next[1];

		const double cost_max_sq = cost_prev_sq >cost_next_sq ? cost_prev_sq : cost_next_sq;
		k_split->heap_node = HEAP_insert(heap, cost_max_sq, c);
	}
	else {
		if (k_split->heap_node) {
			struct KnotCornerState *c;
			c = HEAP_node_ptr(k_split->heap_node);
			HEAP_remove(heap, k_split->heap_node);
			free(c);
			k_split->heap_node = NULL;
		}
	}
}

/**
 * Attempt to collapse close knots into corners,
 * as long as they fall below the error threshold.
 */
static uint curve_incremental_simplify_corners(
        const double *points, const uint points_len,
        struct Knot *knots, const uint knots_len, uint knots_len_remaining,
        const double error_sq_max, const double error_sq_2x_max,
        const double corner_angle,
        const uint dims,
        uint *r_corner_index_len)
{
	Heap *heap = HEAP_new(0);
	double plane_no[dims];

	double k_proj_ref[dims];
	double k_proj_split[dims];

	const double corner_angle_cos = cos(corner_angle);

	uint corner_index_len = 0;

	for (uint i = 0; i < knots_len; i++) {
		if ((knots[i].is_removed == false) &&
		    (knots[i].can_remove == true) &&
		    (knots[i].next && knots[i].next->can_remove))
		{
			struct Knot *k_prev = &knots[i];
			struct Knot *k_next = k_prev->next;

			/* angle outside threshold */
			if (dot_vnvn(k_prev->tan[0], k_next->tan[1], dims) < corner_angle_cos) {
				/* Measure distance projected onto a plane,
				 * since the points may be offset along their own tangents. */
				sub_vn_vnvn(plane_no, k_next->tan[0], k_prev->tan[1], dims);

				/* compare 2x so as to allow both to be changed by maximum of error_sq_max */
				const unsigned int split_knot_index = knot_find_split_point(
				        points, points_len,
				        k_prev, k_next,
				        knots_len,
				        plane_no,
				        dims);

				if (split_knot_index != (unsigned int)-1) {

					project_vn_vnvn(k_proj_ref,   &points[k_prev->point_index * dims], k_prev->tan[1], dims);
					project_vn_vnvn(k_proj_split, &points[split_knot_index    * dims], k_prev->tan[1], dims);

					if (len_squared_vnvn(k_proj_ref, k_proj_split, dims) < error_sq_2x_max) {

						project_vn_vnvn(k_proj_ref,   &points[k_next->point_index * dims], k_next->tan[0], dims);
						project_vn_vnvn(k_proj_split, &points[split_knot_index    * dims], k_next->tan[0], dims);

						if (len_squared_vnvn(k_proj_ref, k_proj_split, dims) < error_sq_2x_max) {

							struct Knot *k_split = &knots[split_knot_index];

							knot_corner_error_recalculate(
							        heap, points, points_len,
							        k_split,
							        k_prev, k_next,
							        error_sq_max,
							        dims);
						}
					}
				}
			}
		}
	}

	while (HEAP_is_empty(heap) == false) {
		struct KnotCornerState *c = HEAP_popmin(heap);

		struct Knot *k_split = &knots[c->knot_index];

		/* remove while collapsing */
		struct Knot *k_prev_remove  = &knots[c->knot_index_adjacent[0]];
		struct Knot *k_next_remove  = &knots[c->knot_index_adjacent[1]];

		struct Knot *k_prev = k_prev_remove->prev;
		struct Knot *k_next = k_next_remove->next;

		if (k_prev == NULL || k_next == NULL) {
			free(c);
			continue;
		}

		/* insert */
		k_split->is_removed = false;
		k_split->prev = k_prev;
		k_split->next = k_next;
		k_prev->next = k_split;
		k_next->prev = k_split;

		/* remove */
		k_prev_remove->is_removed = true;
		k_next_remove->is_removed = true;
		k_prev_remove->prev = NULL;
		k_prev_remove->next = NULL;
		k_next_remove->prev = NULL;
		k_next_remove->next = NULL;

		/* update tangents */
		k_split->tan[0] = k_prev_remove->tan[0];
		k_split->tan[1] = k_next_remove->tan[1];

		/* own handles */
		k_prev->handles[1]  = c->handles_prev[0];
		k_split->handles[0] = c->handles_prev[1];
		k_split->handles[1] = c->handles_next[0];
		k_next->handles[0]  = c->handles_next[1];

		k_split->heap_node = NULL;

		free(c);

#define UPDATE_KNOT(k_update) \
		if ((k_update) && (k_update)->heap_node) { \
			c = HEAP_node_ptr((k_update)->heap_node); \
			knot_corner_error_recalculate( \
			        heap, points, points_len, \
			        (k_update), \
			        &knots[c->knot_index_adjacent[0]], &knots[c->knot_index_adjacent[0]], \
			        error_sq_max, \
			        dims); \
		} ((void)0)

		/* chances are _very_ low these will be corners, its only a NULL check though. */
		UPDATE_KNOT(k_prev);
		UPDATE_KNOT(k_next);

		/* possible these are corners. */
		UPDATE_KNOT(k_prev->prev);
		UPDATE_KNOT(k_next->next);

		// printf("Reducing!\n");

		k_split->is_corner = true;

		knots_len_remaining--;
		corner_index_len++;
	}

	HEAP_free(heap, free);

	*r_corner_index_len = corner_index_len;

	return knots_len_remaining;
}

#endif  /* USE_CORNER_DETECT */


int curve_fit_cubic_to_points_incremental_db(
        const double *points,
        const uint    points_len,
        const uint    dims,
        const double  error_threshold,
        const uint    calc_flag,
        const uint   *corners,
        const uint    corners_len,
        const double  corner_angle,

        double **r_cubic_array, uint *r_cubic_array_len,
        uint   **r_cubic_orig_index,
        uint   **r_corner_index_array, uint *r_corner_index_len)
{
//	const double error_threshold_sq = error_threshold * error_threshold;
	const uint knots_len = points_len;
	Knot *knots = malloc(sizeof(Knot) * knots_len);
	knots[0].next = NULL;

(void)calc_flag;
(void)corners;
(void)corners_len;

	const bool is_cyclic = false;
	const bool use_corner = (corner_angle < M_PI);

	double *tangents = malloc(sizeof(double) * knots_len * 2 * dims);

	{
		double *t_step = tangents;
		for (uint i = 0; i < knots_len; i++) {
			knots[i].next = (knots + i) + 1;
			knots[i].prev = (knots + i) - 1;

			knots[i].heap_node = NULL;
			knots[i].knot_index = i;
			knots[i].point_index = i;
			knots[i].can_remove = true;
			knots[i].is_removed = false;
			knots[i].is_corner = false;
			knots[i].tan[0] = t_step; t_step += dims;
			knots[i].tan[1] = t_step; t_step += dims;
		}
		assert(t_step == &tangents[knots_len * 2 * dims]);
	}

	if (is_cyclic) {
		knots[0].prev = &knots[knots_len - 1];
		knots[knots_len - 1].next = &knots[0];
	}
	else {
		knots[0].prev = NULL;
		knots[knots_len - 1].next = NULL;

		/* always keep end-points */
		knots[0].can_remove = false;
		knots[knots_len - 1].can_remove = false;
	}

	for (uint i = 0; i < knots_len; i++) {
		Knot *k = &knots[i];
		double a[dims];
		double b[dims];

		if (k->prev) {
			sub_vn_vnvn(a, &points[k->prev->point_index * dims], &points[k->point_index * dims], dims);
			normalize_vn(a, dims);
		}
		else {
			zero_vn(a, dims);
		}

		if (k->next) {
			sub_vn_vnvn(b, &points[k->point_index * dims], &points[k->next->point_index * dims], dims);
			normalize_vn(b, dims);
		}
		else {
			zero_vn(b, dims);
		}

		add_vn_vnvn(k->tan[0], a, b, dims);
		normalize_vn(k->tan[0], dims);
		copy_vnvn(k->tan[1], k->tan[0], dims);
	}

	uint knots_len_remaining = knots_len;

	knots_len_remaining = curve_incremental_simplify(
	        points, points_len,
	        knots, knots_len, knots_len_remaining,
	        /* XXX, 'use_corner' works but is weak*/
	        SQUARE(error_threshold / (use_corner ? 2 : 1)), dims);

#ifdef USE_CORNER_DETECT
	if (use_corner) {
		for (uint i = 0; i < knots_len; i++) {
			assert(knots[i].heap_node == NULL);
		}

		knots_len_remaining = curve_incremental_simplify_corners(
		        points, points_len,
		        knots, knots_len, knots_len_remaining,
		        SQUARE(error_threshold), SQUARE(error_threshold * 2),
		        corner_angle,
		        dims,
		        r_corner_index_len);

		/* XXX, needed because of use_corner hack above */
		knots_len_remaining = curve_incremental_simplify(
		        points, points_len,
		        knots, knots_len, knots_len_remaining,
		        SQUARE(error_threshold), dims);

		if (is_cyclic == false) {
			*r_corner_index_len += 2;
		}

		uint *corner_index_array = malloc(sizeof(uint) * (*r_corner_index_len));
		uint k_index = 0, c_index = 0;
		uint i = 0;

		if (is_cyclic == false) {
			corner_index_array[c_index++] = k_index;
			k_index++;
			i++;
		}

		for (; i < knots_len; i++) {
			if (knots[i].is_removed == false) {
				if (knots[i].is_corner == true) {
					corner_index_array[c_index++] = k_index;
				}
				k_index++;
			}
		}

		if (is_cyclic == false) {
			corner_index_array[c_index++] = k_index;
			k_index++;
		}

		assert(c_index == *r_corner_index_len);
		*r_corner_index_array = corner_index_array;
	}
#endif

	uint *cubic_orig_index = NULL;

	if (r_cubic_orig_index) {
		cubic_orig_index = malloc(sizeof(uint) * knots_len_remaining);
	}

	Knot *knots_first = NULL;
	{
		Knot *k;
		for (uint i = 0; i < knots_len; i++) {
			if (knots[i].is_removed == false) {
				knots_first = &knots[i];
				break;
			}
		}

		if (cubic_orig_index) {
			k = knots_first;
			for (uint i = 0; i < knots_len_remaining; i++, k = k->next) {
				cubic_orig_index[i] = k->point_index;
			}
		}
	}

	/* 3x for one knot and two handles */
	double *cubic_array = malloc(sizeof(double) * knots_len_remaining * 3 * dims);

	{
		double *c_step = cubic_array;
		Knot *k = knots_first;
		for (uint i = 0; i < knots_len_remaining; i++, k = k->next) {
			const double *p = &points[k->point_index * dims];

			madd_vn_vnvn_fl(c_step, p, k->tan[0], k->handles[0], dims);
			c_step += dims;
			copy_vnvn(c_step, p, dims);
			c_step += dims;
			madd_vn_vnvn_fl(c_step, p, k->tan[1], k->handles[1], dims);
			c_step += dims;
		}
		assert(c_step == &cubic_array[knots_len_remaining * 3 * dims]);
	}

	free(knots);
	free(tangents);

	if (r_cubic_orig_index) {
		*r_cubic_orig_index = cubic_orig_index;
	}

	*r_cubic_array = cubic_array;
	*r_cubic_array_len = knots_len_remaining;

	return 0;
}


int curve_fit_cubic_to_points_incremental_fl(
        const float          *points,
        const unsigned int    points_len,
        const unsigned int    dims,
        const float           error_threshold,
        const unsigned int    calc_flag,
        const unsigned int   *corners,
        unsigned int          corners_len,
        const float           corner_angle,

        float **r_cubic_array, unsigned int *r_cubic_array_len,
        unsigned int   **r_cubic_orig_index,
        unsigned int   **r_corner_index_array, unsigned int *r_corner_index_len)
{
	const uint points_flat_len = points_len * dims;
	double *points_db = malloc(sizeof(double) * points_flat_len);

	copy_vndb_vnfl(points_db, points, points_flat_len);

	double *cubic_array_db = NULL;
	float  *cubic_array_fl = NULL;
	uint    cubic_array_len = 0;

	int result = curve_fit_cubic_to_points_incremental_db(
	        points_db, points_len, dims, error_threshold, calc_flag, corners, corners_len,
	        corner_angle,
	        &cubic_array_db, &cubic_array_len,
	        r_cubic_orig_index,
	        r_corner_index_array, r_corner_index_len);
	free(points_db);

	if (!result) {
		uint cubic_array_flat_len = cubic_array_len * 3 * dims;
		cubic_array_fl = malloc(sizeof(float) * cubic_array_flat_len);
		for (uint i = 0; i < cubic_array_flat_len; i++) {
			cubic_array_fl[i] = (float)cubic_array_db[i];
		}
		free(cubic_array_db);
	}

	*r_cubic_array = cubic_array_fl;
	*r_cubic_array_len = cubic_array_len;

	return result;
}