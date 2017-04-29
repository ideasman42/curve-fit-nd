
"""
Example of watching a single test:
  watch -n2 "USE_SVG=1 nice -n 20 python -m unittest tests.FreehandTest.test_curve_01"
"""

import os
import sys
import math

# ----------------------------------------------------------------------------
# Math functions


def sub_vnvn(v0, v1):
    return tuple(a - b for a, b in zip(v0, v1))


def dot_vnvn(v0, v1):
    return sum(a * b for a, b in zip(v0, v1))


def len_squared_vn(v0):
    return dot_vnvn(v0, v0)


def len_squared_vnvn(v0, v1):
    d = sub_vnvn(v0, v1)
    return len_squared_vn(d)


def len_vnvn(v0, v1):
    return math.sqrt(len_squared_vnvn(v0, v1))


def interp_vnvn(v0, v1, t):
    s = 1.0 - t
    return tuple((s * a) + (t * b) for a, b in zip(v0, v1))


def interp_cubic_vn(v0, v1, v2, v3, u):
    q0 = interp_vnvn(v0, v1, u)
    q1 = interp_vnvn(v1, v2, u)
    q2 = interp_vnvn(v2, v3, u)

    r0 = interp_vnvn(q0, q1, u)
    r1 = interp_vnvn(q1, q2, u)

    return interp_vnvn(r0, r1, u)


# ----------------------------------------------------------------------------

# refine error calculation
USE_REFINE = True
REFINE_STEPS = 6
REFINE_SHRINK = 0.5

import unittest

# module from ../c_python_ext
import curve_fit_nd

TEST_DATA_PATH = os.path.join(os.path.dirname(__file__), "data")

USE_SVG = os.environ.get("USE_SVG")

sys.path.append(TEST_DATA_PATH)


def iter_pairs(iterable):
    it = iter(iterable)
    prev_ = next(it)
    while True:
        try:
            next_ = next(it)
        except StopIteration:
            return
        yield (prev_, next_)
        prev_ = next_


def export_svg(name, s, points, measure_points):
    # write svg's to tests dir for now
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
           'width="%d" height="%d" '
           'viewBox="%d %d %d %d" '
           'xmlns="http://www.w3.org/2000/svg">\n' %
           (
            # width, height
            int(margin * scale * 2), int(margin * scale * 2),
            # viewBox
            -int(margin * scale), -int(margin * scale), int(margin * scale * 2), int(margin * scale * 2),
            ))

        fw('<rect x="%d" y="%d" width="%d" height="%d" fill="black"/>\n' %
            (-int(margin * scale), -int(margin * scale), int(margin * scale * 2), int(margin * scale * 2),
            ))

        if s:
            fw('<g stroke="white" stroke-opacity="0.25" stroke-width="1">\n')
            for (i0, v0), (i1, v1) in iter_pairs(s):
                k0 = v0[1][0] * scale, -v0[1][1] * scale
                h0 = v0[2][0] * scale, -v0[2][1] * scale

                h1 = v1[0][0] * scale, -v1[0][1] * scale
                k1 = v1[1][0] * scale, -v1[1][1] * scale

                fw('<path d="M%.4f,%.4f C%.4f,%.4f %.4f,%.4f %.4f,%.4f" />\n' %
                   (k0[0], k0[1],
                    h0[0], h0[1],
                    h1[0], h1[1],
                    k1[0], k1[1],
                    ))
            fw('</g>\n')
        if points:
            fw('<g fill="yellow" fill-opacity="0.25" stroke="white" stroke-opacity="0.25" stroke-width="1">\n')
            for v0 in points:
                fw('<circle cx="%.4f" cy="%.4f" r="1"/>\n' %
                (v0[0] * scale, -v0[1] * scale))
            fw('</g>\n')

        if s:
            # tangent handles
            fw('<g fill="white" fill-opacity="0.5" stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for i, p in s:
                for v in p:
                    fw('<circle cx="%.4f" cy="%.4f" r="2"/>\n' %
                    (v[0] * scale, -v[1] * scale))
            fw('</g>\n')

            fw('<g fill="white" fill-opacity="0.5" stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for i, p in s:
                v = p[1]
                fw('<circle cx="%.4f" cy="%.4f" r="2"/>\n' %
                (v[0] * scale, -v[1] * scale))
            fw('</g>\n')


            # lines
            fw('<g stroke="white" stroke-opacity="0.2" stroke-width="1">\n')
            for i, (v0, v1, v2) in s:
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v0[0] * scale, -v0[1] * scale, v1[0] * scale, -v1[1] * scale))
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v1[0] * scale, -v1[1] * scale, v2[0] * scale, -v2[1] * scale))
            fw('</g>\n')

        if measure_points:
            fw('<g stroke="white" stroke-opacity="0.5" stroke-width="0.5">\n')
            for e0, e1 in measure_points:
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (e0[0] * scale, -e0[1] * scale, e1[0] * scale, -e1[1] * scale))
            fw('</g>\n')

        fw('</svg>')


def curve_fit(s, error, corner_angle, is_cyclic):
    if corner_angle is None:
        corner_angle = 10
    c = curve_fit_nd.curve_from_points(s, error, corner_angle, is_cyclic)
    return c


def curve_error_max(points, curve, r_measure_points):
    error_max_sq = 0.0
    for (i0, v0), (i1, v1) in iter_pairs(curve):
        k0 = v0[1]
        h0 = v0[2]

        h1 = v1[0]
        k1 = v1[1]

        lens = [None] * (i1 - i0)
        lens_accum = 0.0
        for i, (s0, s1) in enumerate(iter_pairs(range(i0, i1 + 1))):
            lens[i] = lens_accum
            lens_accum += len_vnvn(points[s0], points[s1])

        do_refine = USE_REFINE
        if lens[-1] != 0.0:
            lens[:] = [l / lens[-1] for l in lens]
        else:
            do_refine = False

        for i, s in enumerate(range(i0, i1)):
            u = lens[i]
            p_real = points[s]
            p_curve = interp_cubic_vn(k0, h0, h1, k1, u)

            # step up and down the cubic to reach a close point
            if do_refine:
                error_best = len_squared_vnvn(p_real, p_curve)

                u_step = 0.0
                u_tot = 0.0
                if i != 0:
                    u_step += lens[i - 1] - u
                    # assert(u - lens[i - 1] >= 0.0)
                    u_tot += 1.0
                if i != len(lens) - 1:
                    u_step += u - lens[i + 1]
                    # assert(lens[i + 1] - u >= 0.0)
                    u_tot += 1.0

                u_step /= u_tot

                # refine, start at half
                u_step *= REFINE_SHRINK

                u_best = u
                p_best = p_curve
                refine_count = 0


                while True:
                    u_init = u_best
                    for u_test in (u_best + u_step, u_best - u_step):
                        p_test = interp_cubic_vn(k0, h0, h1, k1, u_test)
                        error_test = len_squared_vnvn(p_real, p_test)
                        if error_test < error_best:
                            error_best = error_test
                            u_best = u_test
                            p_best = p_test

                    if u_init == u_best:
                        refine_count += 1
                        if refine_count == REFINE_STEPS:
                            break
                        else:
                            # refine further
                            u_step *= REFINE_SHRINK
                error_max_sq = max(error_max_sq, error_best)
                p_curve = p_best
                # end USE_REFINE

            else:
                error_max_sq = max(error_max_sq, len_squared_vnvn(p_real, p_curve))


            r_measure_points.append((p_real, p_curve))

    return math.sqrt(error_max_sq)


def test_data_load(name):
    return __import__(name).data


class TestDataFile_Helper:

    def assertTestData(self, name, error, corner_angle=None, is_cyclic=False):
        points = test_data_load(name)

        curve = curve_fit(points, error, corner_angle, is_cyclic)
        # print(name + ix_id)

        measure_points = []
        error_test = curve_error_max(points, curve, measure_points)
        # print(error_test, error)

        if USE_SVG:
            export_svg(name, curve, points, measure_points)

        # scale the error up to allow for some minor discrepancy in USE_REFINE
        self.assertLess(error_test, error * 1.01)


class FreehandTest(unittest.TestCase, TestDataFile_Helper):
    """
    Test freehand curve fitting.
    """

    def test_curve_01(self):
        self.assertTestData("test_curve_freehand_01", 0.01)

    def test_curve_02(self):
        self.assertTestData("test_curve_freehand_02", 0.01)

    def test_curve_03(self):
        self.assertTestData("test_curve_freehand_03", 0.01, math.radians(30))

    def test_curve_04(self):
        self.assertTestData("test_curve_freehand_04_cyclic", 0.0075, math.radians(70), is_cyclic=True)


if __name__ == "__main__":
    unittest.main()
