"""
Example of watching a single test:
  watch -n2 "USE_SVG=1 nice -n 20 python -m unittest tests.CurveFitTest.test_curve_01"
"""

import glob
import math
import os
import sys
import unittest

from importlib import import_module
from itertools import pairwise
from typing import Callable, NamedTuple, Sequence, cast

import curve_fit_nd  # type: ignore[import-not-found]

# Type aliases for points/vectors.
float2 = tuple[float, float]
floatN = tuple[float, ...]

# ----------------------------------------------------------------------------
# Math functions


def sub_vnvn(v0: Sequence[float], v1: Sequence[float]) -> floatN:
    return tuple(a - b for a, b in zip(v0, v1))


def dot_vnvn(v0: Sequence[float], v1: Sequence[float]) -> float:
    return sum(a * b for a, b in zip(v0, v1))


def len_squared_vn(v0: Sequence[float]) -> float:
    return dot_vnvn(v0, v0)


def len_squared_vnvn(v0: Sequence[float], v1: Sequence[float]) -> float:
    d = sub_vnvn(v0, v1)
    return len_squared_vn(d)


def len_vnvn(v0: Sequence[float], v1: Sequence[float]) -> float:
    return math.sqrt(len_squared_vnvn(v0, v1))


def cross_v2v2(v0: float2, v1: float2) -> float:
    """2D cross product (z-component of 3D cross product)."""
    return v0[0] * v1[1] - v0[1] * v1[0]


def interp_vnvn(v0: Sequence[float], v1: Sequence[float], t: float) -> floatN:
    s = 1.0 - t
    return tuple((s * a) + (t * b) for a, b in zip(v0, v1))


def interp_cubic_vn(
        v0: Sequence[float],
        v1: Sequence[float],
        v2: Sequence[float],
        v3: Sequence[float],
        u: float,
) -> floatN:
    q0 = interp_vnvn(v0, v1, u)
    q1 = interp_vnvn(v1, v2, u)
    q2 = interp_vnvn(v2, v3, u)

    r0 = interp_vnvn(q0, q1, u)
    r1 = interp_vnvn(q1, q2, u)

    return interp_vnvn(r0, r1, u)


def reflect_vnvn(v0: floatN, v1: floatN) -> floatN:
    """Reflect v0 across v1."""
    return tuple(2.0 * b - a for a, b in zip(v0, v1))


def interp_catmull_rom_midpoint_vnvnvnvn(
        p0: floatN,
        p1: floatN,
        p2: floatN,
        p3: floatN,
) -> floatN:
    """
    Compute smooth midpoint between p1 and p2 using Catmull-Rom interpolation.

    Uses the 4-point subdivision rule: (-p0 + 9*p1 + 9*p2 - p3) / 16.
    """
    return tuple((-a + 9.0 * b + 9.0 * c - d) / 16.0 for a, b, c, d in zip(p0, p1, p2, p3))


def subdivide_points(
        points: Sequence[floatN],
        iterations: int = 1,
        is_cyclic: bool = False,
) -> list[floatN]:
    """
    Subdivide points by inserting smooth midpoints between each pair.

    Uses 4-point Catmull-Rom subdivision for smooth curvature-following placement.
    Original points are preserved at their exact positions.

    :arg points: Input points (n-dimensional).
    :arg iterations: Number of subdivision passes (each roughly doubles point count).
    :arg is_cyclic: If True, also subdivide the segment from last to first point.
    :return: Subdivided points with original points preserved.
    """
    result: list[floatN] = list(points)

    for _ in range(iterations):
        n = len(result)
        if n < 2:
            break

        num_segments = n if is_cyclic else n - 1
        # Pre-allocate: each segment adds 1 midpoint, plus original points.
        # Cyclic: n original + n midpoints = 2n.
        # Non-cyclic: n original + (n-1) midpoints = 2n - 1.
        new_len = 2 * n if is_cyclic else 2 * n - 1
        new_points: list[floatN] = [None] * new_len  # type: ignore[list-item]

        for i in range(num_segments):
            # Place original point at even index.
            new_points[i * 2] = result[i]

            # Get 4 points for interpolation.
            if is_cyclic:
                p0 = result[(i - 1) % n]
                p1 = result[i]
                p2 = result[(i + 1) % n]
                p3 = result[(i + 2) % n]
            else:
                p1 = result[i]
                p2 = result[i + 1]
                # Reflect for boundary cases.
                p0 = result[i - 1] if i > 0 else reflect_vnvn(p2, p1)
                p3 = result[i + 2] if i + 2 < n else reflect_vnvn(p1, p2)

            # Place midpoint at odd index.
            new_points[i * 2 + 1] = interp_catmull_rom_midpoint_vnvnvnvn(p0, p1, p2, p3)

        # Add final point for non-cyclic.
        if not is_cyclic:
            new_points[-1] = result[-1]

        result = new_points

    return result


# ----------------------------------------------------------------------------
# Constants

# Minimum subdivisions per Bezier segment for area calculation.
AREA_MIN_SUBDIVISIONS = 16

# Refinement settings for finding closest curve point.
USE_REFINE = True
REFINE_STEPS = 8
REFINE_SHRINK = 0.5

# Tolerance scale for error comparison (allows minor discrepancy in refinement).
ERROR_TOLERANCE_SCALE = 1.01

# Tolerance for area delta comparison (accounts for cross-platform float differences).
AREA_DELTA_TOLERANCE = 1e-9

# Default corner angle in radians (about 10 degrees).
DEFAULT_CORNER_ANGLE = math.radians(10)

TEST_DATA_PATH = os.path.join(os.path.dirname(__file__), "data")

USE_SVG = os.environ.get("USE_SVG") or False
USE_STRESS_TEST = os.environ.get("USE_STRESS_TEST") or False

# Number of subdivision iterations for stress testing.
STRESS_TEST_SUBDIV_ITERATIONS = 5

if USE_SVG:
    svg_dir = os.path.join(TEST_DATA_PATH, "..", "data_svg")
    if os.path.isdir(svg_dir):
        for svg_file in glob.glob(os.path.join(svg_dir, "*.svg")):
            os.remove(svg_file)

sys.path.append(TEST_DATA_PATH)


def export_svg(
        name: str,
        s: Sequence[tuple[int, tuple[float2, float2, float2]]],
        points: Sequence[float2],
        measure_points: Sequence[tuple[float2, float2]],
) -> None:
    # Write SVGs to tests dir for now.
    dirname = os.path.join(TEST_DATA_PATH, "..", "data_svg")
    os.makedirs(dirname, exist_ok=True)
    fp = os.path.join(dirname, name + ".svg")
    scale = 1024.0
    margin = 1.1
    with open(fp, 'w', encoding="utf-8") as f:
        fw = f.write
        fw('<?xml version="1.0" encoding="UTF-8"?>\n')
        fw('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1 Tiny//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11-tiny.dtd">\n')
        fw('<svg version="1.1" '
           'width="{:d}" height="{:d}" '
           'viewBox="{:d} {:d} {:d} {:d}" '
           'xmlns="http://www.w3.org/2000/svg">\n'.format(
               # width, height
               int(margin * scale * 2), int(margin * scale * 2),
               # viewBox
               -int(margin * scale), -int(margin * scale), int(margin * scale * 2), int(margin * scale * 2),
           ))

        fw('<rect x="{:d}" y="{:d}" width="{:d}" height="{:d}" fill="black"/>\n'.format(
            -int(margin * scale), -int(margin * scale), int(margin * scale * 2), int(margin * scale * 2),
        ))

        if s:
            fw('<g stroke="white" stroke-opacity="0.25" stroke-width="1">\n')
            for (_i0, v0), (_i1, v1) in pairwise(s):
                k0 = v0[1][0] * scale, -v0[1][1] * scale
                h0 = v0[2][0] * scale, -v0[2][1] * scale

                h1 = v1[0][0] * scale, -v1[0][1] * scale
                k1 = v1[1][0] * scale, -v1[1][1] * scale

                fw('<path d="M{:.4f},{:.4f} C{:.4f},{:.4f} {:.4f},{:.4f} {:.4f},{:.4f}" />\n'.format(
                    k0[0], k0[1],
                    h0[0], h0[1],
                    h1[0], h1[1],
                    k1[0], k1[1],
                ))
            fw('</g>\n')
        if points:
            fw('<g fill="yellow" fill-opacity="0.25" stroke="white" stroke-opacity="0.25" stroke-width="1">\n')
            for pt in points:
                fw('<circle cx="{:.4f}" cy="{:.4f}" r="1"/>\n'.format(
                    pt[0] * scale, -pt[1] * scale))
            fw('</g>\n')

        if s:
            # Tangent handles.
            fw('<g fill="white" fill-opacity="0.5" stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for _i, (hl, _co, hr) in s:
                for v in (hl, hr):
                    fw('<circle cx="{:.4f}" cy="{:.4f}" r="2"/>\n'.format(
                        v[0] * scale, -v[1] * scale))
            fw('</g>\n')

            # Knots.
            fw('<g fill="white" fill-opacity="0.5" stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for _i, (_hl, co, _hr) in s:
                fw('<circle cx="{:.4f}" cy="{:.4f}" r="2"/>\n'.format(
                    co[0] * scale, -co[1] * scale))
            fw('</g>\n')

            # Lines.
            fw('<g stroke="white" stroke-opacity="0.2" stroke-width="1">\n')
            for _i, (hl, co, hr) in s:
                fw('<line x1="{:.4f}" y1="{:.4f}" x2="{:.4f}" y2="{:.4f}" />\n'.format(
                    hl[0] * scale, -hl[1] * scale, co[0] * scale, -co[1] * scale))
                fw('<line x1="{:.4f}" y1="{:.4f}" x2="{:.4f}" y2="{:.4f}" />\n'.format(
                    co[0] * scale, -co[1] * scale, hr[0] * scale, -hr[1] * scale))
            fw('</g>\n')

        if measure_points:
            fw('<g stroke="white" stroke-opacity="0.5" stroke-width="0.5">\n')
            for e0, e1 in measure_points:
                fw('<line x1="{:.4f}" y1="{:.4f}" x2="{:.4f}" y2="{:.4f}" />\n'.format(
                    e0[0] * scale, -e0[1] * scale, e1[0] * scale, -e1[1] * scale))
            fw('</g>\n')

        fw('</svg>')


def curve_fit(
        points: Sequence[floatN],
        error: float,
        corner_angle: float | None,
        is_cyclic: bool,
) -> list[tuple[int, tuple[float2, float2, float2]]]:
    if corner_angle is None:
        corner_angle = DEFAULT_CORNER_ANGLE
    return curve_fit_nd.curve_from_points(points, error, corner_angle, is_cyclic)  # type: ignore[no-any-return]


def interp_polyline_at_arc_length(
        points: Sequence[float2],
        indices: list[int],
        arc_lens: list[float],
        u: float,
) -> float2:
    """
    Interpolate along a polyline at normalized arc-length parameter u.

    :arg points: Full points array.
    :arg indices: Indices into points forming the polyline segment.
    :arg arc_lens: Normalized arc-lengths for each index (0.0 to 1.0).
    :arg u: Parameter in range [0, 1].
    :return: Interpolated point on the polyline.
    """
    if u <= 0.0:
        return points[indices[0]]
    if u >= 1.0:
        return points[indices[-1]]

    # Find which segment contains u
    for k in range(len(arc_lens) - 1):
        if arc_lens[k] <= u <= arc_lens[k + 1]:
            # Interpolate within this segment
            segment_start = arc_lens[k]
            segment_end = arc_lens[k + 1]
            if segment_end == segment_start:
                t = 0.0
            else:
                t = (u - segment_start) / (segment_end - segment_start)
            return cast(float2, interp_vnvn(points[indices[k]], points[indices[k + 1]], t))

    # Fallback (shouldn't reach here)
    return points[indices[-1]]


def curve_aligned_samples(
        points: Sequence[float2],
        curve: Sequence[tuple[int, tuple[float2, float2, float2]]],
        is_cyclic: bool,
        min_subdivisions: int,
        refine_threshold: float,
) -> list[tuple[float2, float2]]:
    """
    Generate aligned (original_point, curve_point) pairs.

    Each segment gets at least min_subdivisions samples.
    Adaptively adds more where deviation exceeds refine_threshold.

    :arg points: Original input points.
    :arg curve: Fitted curve data.
    :arg is_cyclic: Whether the curve is cyclic.
    :arg min_subdivisions: Minimum samples per segment.
    :arg refine_threshold: Max deviation before adding samples (0 disables).
    :return: List of (original_point, curve_point) tuples.
    """
    n = len(points)
    result: list[tuple[float2, float2]] = []

    # Build segment list (with wrap-around for cyclic)
    segments = list(pairwise(curve))
    if is_cyclic:
        segments.append((curve[-1], curve[0]))

    for (i0, v0), (i1, v1) in segments:
        k0, h0 = v0[1], v0[2]
        h1, k1 = v1[0], v1[1]

        # Handle index wrap-around for cyclic curves
        if i1 > i0:
            indices = list(range(i0, i1 + 1))
        else:
            indices = list(range(i0, n)) + list(range(0, i1 + 1))

        # Compute arc length parameterization for original polyline
        arc_lens = [0.0]
        for index0, index1 in pairwise(indices):
            arc_lens.append(arc_lens[-1] + len_vnvn(points[index0], points[index1]))

        total_length = arc_lens[-1]
        if total_length > 0.0:
            arc_lens = [length / total_length for length in arc_lens]

        def sample_at(u: float) -> tuple[float2, float2]:
            """Sample both curves at parameter u."""
            curve_pt = cast(float2, interp_cubic_vn(k0, h0, h1, k1, u))
            orig_pt = interp_polyline_at_arc_length(points, indices, arc_lens, u)
            return (orig_pt, curve_pt)

        def deviation(p0: tuple[float2, float2], p1: tuple[float2, float2], p_mid: tuple[float2, float2]) -> float:
            """Measure how far midpoint deviates from linear interpolation."""
            # Interpolate expected midpoint for both orig and curve
            orig_expected = interp_vnvn(p0[0], p1[0], 0.5)
            curve_expected = interp_vnvn(p0[1], p1[1], 0.5)
            # Measure deviation as max of both
            orig_dev = len_vnvn(orig_expected, p_mid[0])
            curve_dev = len_vnvn(curve_expected, p_mid[1])
            return max(orig_dev, curve_dev)

        def subdivide(u_start: float, u_end: float, p_start: tuple[float2, float2], p_end: tuple[float2, float2], depth: int) -> list[tuple[float, tuple[float2, float2]]]:
            """Recursively subdivide if deviation exceeds threshold."""
            u_mid = (u_start + u_end) / 2.0
            p_mid = sample_at(u_mid)

            # Check if we need to subdivide further
            if refine_threshold > 0.0 and depth < 10:  # Max depth to prevent infinite recursion
                dev = deviation(p_start, p_end, p_mid)
                if dev > refine_threshold:
                    # Subdivide both halves
                    left = subdivide(u_start, u_mid, p_start, p_mid, depth + 1)
                    right = subdivide(u_mid, u_end, p_mid, p_end, depth + 1)
                    return left + [(u_mid, p_mid)] + right

            return [(u_mid, p_mid)]

        # Start with uniform samples
        segment_samples: list[tuple[float, tuple[float2, float2]]] = []
        for sub in range(min_subdivisions + 1):
            u = sub / min_subdivisions
            segment_samples.append((u, sample_at(u)))

        # Adaptively refine between initial samples
        if refine_threshold > 0.0:
            refined: list[tuple[float, tuple[float2, float2]]] = []
            for (u0, p0), (u1, p1) in pairwise(segment_samples):
                refined.append((u0, p0))
                refined.extend(subdivide(u0, u1, p0, p1, 0))
            refined.append(segment_samples[-1])
            segment_samples = refined

        # Sort by parameter (only needed when refinement adds interleaved points)
        if refine_threshold > 0.0:
            segment_samples.sort(key=lambda x: x[0])

        # Add to result (skip first point if not first segment to avoid duplicates)
        start_idx = 0 if not result else 1
        for _, sample in segment_samples[start_idx:]:
            result.append(sample)

    return result


def curve_error_max(
        points: Sequence[float2],
        curve: Sequence[tuple[int, tuple[float2, float2, float2]]],
        is_cyclic: bool,
) -> tuple[float, list[tuple[float2, float2]]]:
    """
    Return (max_error, measure_points) for visualization.

    For each original point, finds the closest point on the Bezier curve.
    """
    n = len(points)
    error_max_sq = 0.0
    measure_points: list[tuple[float2, float2]] = []

    # Build segment list (with wrap-around for cyclic)
    segments = list(pairwise(curve))
    if is_cyclic:
        segments.append((curve[-1], curve[0]))

    for (i0, v0), (i1, v1) in segments:
        k0, h0 = v0[1], v0[2]
        h1, k1 = v1[0], v1[1]

        # Handle index wrap-around for cyclic curves
        if i1 > i0:
            indices = list(range(i0, i1))
        else:
            indices = list(range(i0, n)) + list(range(0, i1))

        # Compute arc length parameterization
        arc_lens = [0.0]
        for j in range(len(indices)):
            if j > 0:
                arc_lens.append(arc_lens[-1] + len_vnvn(points[indices[j-1]], points[indices[j]]))

        total_len = arc_lens[-1] if arc_lens else 0.0
        do_refine = USE_REFINE and total_len != 0.0
        if total_len != 0.0:
            arc_lens = [length / total_len for length in arc_lens]

        for j, idx in enumerate(indices):
            u = arc_lens[j] if j < len(arc_lens) else 1.0
            p_real = points[idx]
            p_curve = cast(float2, interp_cubic_vn(k0, h0, h1, k1, u))

            # Refine to find closest point on curve
            if do_refine:
                error_best = len_squared_vnvn(p_real, p_curve)

                # Compute step size from neighboring arc lengths
                u_step = 0.0
                u_tot = 0.0
                if j != 0:
                    u_step += arc_lens[j - 1] - u
                    u_tot += 1.0
                if j != len(arc_lens) - 1:
                    u_step += u - arc_lens[j + 1] if j + 1 < len(arc_lens) else 0.0
                    u_tot += 1.0

                if u_tot > 0.0:
                    u_step /= u_tot
                    u_step *= REFINE_SHRINK

                    u_best = u
                    p_best = p_curve
                    refine_count = 0

                    while True:
                        u_init = u_best
                        for u_test in (u_best + u_step, u_best - u_step):
                            if 0.0 <= u_test <= 1.0:
                                p_test = cast(float2, interp_cubic_vn(k0, h0, h1, k1, u_test))
                                error_test = len_squared_vnvn(p_real, p_test)
                                if error_test < error_best:
                                    error_best = error_test
                                    u_best = u_test
                                    p_best = p_test

                        if u_init == u_best:
                            refine_count += 1
                            if refine_count == REFINE_STEPS:
                                break
                            u_step *= REFINE_SHRINK

                    error_max_sq = max(error_max_sq, error_best)
                    p_curve = p_best
                else:
                    error_max_sq = max(error_max_sq, len_squared_vnvn(p_real, p_curve))
            else:
                error_max_sq = max(error_max_sq, len_squared_vnvn(p_real, p_curve))

            measure_points.append((p_real, p_curve))

    return math.sqrt(error_max_sq), measure_points


def curve_area_delta(
        points: Sequence[float2],
        curve: Sequence[tuple[int, tuple[float2, float2, float2]]],
        is_cyclic: bool,
) -> float:
    """
    Calculate the ribbon area between original points and fitted curve.

    Uses the diagonal cross product method for quadrilateral areas.
    """
    samples = curve_aligned_samples(
        points, curve, is_cyclic,
        min_subdivisions=AREA_MIN_SUBDIVISIONS,
        refine_threshold=0.0,
    )

    area = 0.0
    for (orig_pt, curve_pt), (next_orig_pt, next_curve_pt) in pairwise(samples):
        # Form quadrilateral and compute area
        # Quad: (orig, next_orig, next_curve, curve)
        # Area = 0.5 * |cross(orig - next_curve, next_orig - curve)|
        d1 = cast(float2, sub_vnvn(orig_pt, next_curve_pt))
        d2 = cast(float2, sub_vnvn(next_orig_pt, curve_pt))
        area += abs(cross_v2v2(d1, d2))

    return 0.5 * area


def test_data_load(name: str) -> Sequence[float2]:
    data: Sequence[float2] = import_module(name).data
    return data


def curve_knots_sorted(
        curve: Sequence[tuple[int, tuple[floatN, floatN, floatN]]],
) -> list[tuple[floatN, floatN, floatN]]:
    """
    Extract and rotate curve knots for order-independent comparison.

    Finds the minimum knot by position (co) and rotates the list
    so it comes first, preserving the original order of elements.
    """
    knots = [v for _i, v in curve]
    if not knots:
        return knots
    # Find the index of the minimum knot, ordering by (co, handle_left, handle_right).
    # The knot position (co) is ordered first since it's the most stable identifier,
    # handles depend on neighboring points and may vary with rotation.
    # A stable identifier makes troubleshooting differences between arrays
    # easier since values are more likely to be aligned (though not guaranteed).
    min_index = min(range(len(knots)), key=lambda i: (knots[i][1], knots[i][0], knots[i][2]))
    # Rotate so the minimum is first.
    return knots[min_index:] + knots[:min_index]


def wrap_points_at_middle(points: Sequence[float2]) -> list[float2]:
    """
    Rotate points by half their length.

    For a cyclic curve, this should produce an equivalent curve
    when fitted, proving order-independence.
    """
    n = len(points)
    mid = n // 2
    return list(points[mid:]) + list(points[:mid])


class TestDataFile_MixIn:

    def assertTestData(  # type: ignore[misc]
        self: "CurveFitTest",
        *,
        name: str,
        error: float,
        corner_angle: float | None = None,
        is_cyclic: bool = False,
        expected_knot_count: int = 0,
        expected_area_delta: float = 0.0,
    ) -> None:
        points = test_data_load(name)

        if USE_STRESS_TEST:
            points = cast(
                Sequence[float2],
                subdivide_points(points, STRESS_TEST_SUBDIV_ITERATIONS, is_cyclic),
            )

        curve = curve_fit(points, error, corner_angle, is_cyclic)

        error_test, measure_points = curve_error_max(points, curve, is_cyclic)
        area_delta = curve_area_delta(points, curve, is_cyclic)

        if USE_SVG:
            export_svg(name, curve, points, measure_points)

        self.assertLess(error_test, error * ERROR_TOLERANCE_SCALE)
        if not USE_STRESS_TEST:
            self.assertEqual(len(curve), expected_knot_count)
            self.assertAlmostEqual(area_delta, expected_area_delta, delta=AREA_DELTA_TOLERANCE)

        if is_cyclic:
            self.assertCyclicOrderIndependence(points, curve, error, corner_angle)

    def assertCyclicOrderIndependence(  # type: ignore[misc]
        self: "CurveFitTest",
        points: Sequence[float2],
        curve: Sequence[tuple[int, tuple[float2, float2, float2]]],
        error: float,
        corner_angle: float | None,
    ) -> None:
        """Verify order-independence by wrapping points at middle and comparing."""
        points_wrapped = wrap_points_at_middle(points)
        curve_wrapped = curve_fit(points_wrapped, error, corner_angle, is_cyclic=True)

        # Both should produce the same number of knots.
        self.assertEqual(len(curve_wrapped), len(curve))

        # Sort both curves by knot position and compare.
        knots_orig = curve_knots_sorted(curve)
        knots_wrapped = curve_knots_sorted(curve_wrapped)

        for i, (orig, wrapped) in enumerate(zip(knots_orig, knots_wrapped)):
            # Compare handle_left, co, handle_right.
            for j, (v_orig, v_wrapped) in enumerate(zip(orig, wrapped)):
                self.assertAlmostEqual(
                    v_orig[0], v_wrapped[0],
                    places=9,
                    msg=f"Knot {i} element {j} x mismatch",
                )
                self.assertAlmostEqual(
                    v_orig[1], v_wrapped[1],
                    places=9,
                    msg=f"Knot {i} element {j} y mismatch",
                )


class TestData(NamedTuple):
    # Filename without extension.
    filename: str
    error_max: float
    corner_angle: float | None
    is_cyclic: bool
    expected_knot_count: int
    expected_area_delta: float


test_data = (
    TestData(
        filename="test_curve_freehand_01",
        error_max=0.01,
        corner_angle=None,
        is_cyclic=False,
        expected_knot_count=27,
        expected_area_delta=0.015039482825614363,
    ),
    TestData(
        filename="test_curve_freehand_02",
        error_max=0.01,
        corner_angle=None,
        is_cyclic=False,
        expected_knot_count=30,
        expected_area_delta=0.015650545848830376,
    ),
    TestData(
        filename="test_curve_freehand_03",
        error_max=0.01,
        corner_angle=math.radians(30),
        is_cyclic=False,
        expected_knot_count=20,
        expected_area_delta=0.012887548579986443,
    ),
    TestData(
        filename="test_curve_freehand_04_cyclic",
        error_max=0.0075,
        corner_angle=math.radians(70),
        is_cyclic=True,
        expected_knot_count=28,
        expected_area_delta=0.01025397509596632,
    ),
)


class CurveFitTest(unittest.TestCase, TestDataFile_MixIn):
    pass


def _make_test_from_data(td: TestData) -> Callable[[CurveFitTest], None]:
    def test_method(self: CurveFitTest) -> None:
        self.assertTestData(
            name=td.filename,
            error=td.error_max,
            corner_angle=td.corner_angle,
            is_cyclic=td.is_cyclic,
            expected_knot_count=td.expected_knot_count,
            expected_area_delta=td.expected_area_delta,
        )
    return test_method


for td in test_data:
    setattr(
        CurveFitTest,
        "test_{:s}".format(td.filename),
        _make_test_from_data(td),
    )


if __name__ == "__main__":
    unittest.main()
