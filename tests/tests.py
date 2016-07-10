
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


def closest_to_line_vn(p, l1, l2):
    u = sub_vnvn(l2, l1)
    h = sub_v3vn(p, l1)
    l = dot_vnvn(u, h) / dot_vnvn(u, u)
    cp = tuple(l1_axis + (u_axis * l) for l1_axis, u_axis in zip(l1, u))
    return cp, l


def closest_to_line_segment_vn(p, l1, l2):
    cp, fac = closest_to_line_vn(p, l1, l2)
    # flip checks for !finite case (when segment is a point)
    if not (fac > 0.0):
        return l1
    elif not (fac < 1.0):
        return l2
    else:
        return cp


def dist_squared_to_line_segment_vn(p, l1, l2):
    closest = closest_to_line_segment_vn(closest, p, l1, l2);
    return len_squared_vnvn(closest, p)


def interp_vnvn(v0, v1, t):
    s = 1.0 - t
    return tuple(s * a + t * b for a, b in zip(v0, v1))


def interp_cubic_vn(v0, v1, v2, v3, u):
    q0 = interp_vnvn(v0, v1, u)
    q1 = interp_vnvn(v1, v2, u)
    q2 = interp_vnvn(v2, v3, u)

    r0 = interp_vnvn(q0, q1, u)
    r1 = interp_vnvn(q1, q2, u)

    return interp_vnvn(r0, r1, u)


# ----------------------------------------------------------------------------


import unittest

# module from ../c_python_ext
import curve_fit_nd

TEST_DATA_PATH = os.path.join(os.path.dirname(__file__), "data")

USE_SVG = os.environ.get("USE_SVG")

sys.path.append(TEST_DATA_PATH)


def iter_pairs(iterable):
    '''
    it = iter(iterable)
    prev_ = next(it)
    while True:
        try:
            next_ = next(it)
        except StopIteration:
            return
        yield (prev_, next_)
        prev_ = next_
    '''
    for i in range(1, len(iterable)):
        yield (iterable[i - 1], iterable[i])


def export_svg(name, s, points):
    # write svg's to tests dir for now
    dirname = os.path.join(TEST_DATA_PATH, "..", "data_svg")
    os.makedirs(dirname, exist_ok=True)
    fp = os.path.join(dirname, name + ".svg")
    scale = 512.0
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

                k1 = v1[0][0] * scale, -v1[0][1] * scale
                h1 = v1[1][0] * scale, -v1[1][1] * scale

                fw('<path d="M%.4f,%.4f C%.4f,%.4f %.4f,%.4f %.4f,%.4f" />\n' %
                   (k0[0], k0[1],
                    h0[0], h0[1],
                    k1[0], k1[1],
                    h1[0], h1[1],
                    ))
                # <path d="M100,250 C100,100 400,100 400,250" />
            fw('</g>\n')
        if points:
            fw('<g fill="yellow" fill-opacity="0.25" stroke="white" stroke-opacity="0.25" stroke-width="1">\n')
            for v0 in points:
                fw('<circle cx="%.4f" cy="%.4f" r="1"/>\n' %
                (v0[0] * scale, -v0[1] * scale))
            fw('</g>\n')

        if s:
            fw('<g fill="white" fill-opacity="0.5" stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for i, p in s:
                for v in p:
                    fw('<circle cx="%.4f" cy="%.4f" r="2"/>\n' %
                    (v[0] * scale, -v[1] * scale))

            # lines
            fw('<g stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for i, (v0, v1, v2) in s:
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v0[0] * scale, -v0[1] * scale, v1[0] * scale, -v1[1] * scale))
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v1[0] * scale, -v1[1] * scale, v2[0] * scale, -v2[1] * scale))
            fw('</g>\n')


            fw('</g>\n')


        fw('</svg>')


def dist_squared_to_line_vn(p, points):
    # warning, slow!
    m = len_squared_vnvn(points[0], p)
    for v0, v1 in iter_pairs(points):
        m = min(dist_squared_to_line_segment_vn(p, v1, v2))
    return m


def curve_error_calc(points, curve):

    for v0, v1 in iter_pairs(s):
        k0 = v0[1], -v0[1]
        h0 = v0[2], -v0[2]

        k1 = v1[0], -v1[0]
        h1 = v1[1], -v1[1]

        samples = 8
        for i in range(1, samples):
            fac = i / samples

            v_sample = interp_cubic_vn(k0, h0, h1, k1)



def curve_fit_compare(s, error):
    c = curve_fit_nd.curve_from_points(s, error)
    return c


def test_data_load(name):
    return __import__(name).data


class TestDataFile_Helper:

    def assertTestData(self, name, error):
        points = test_data_load(name)
        curve = curve_fit_compare(points, error)
        if USE_SVG:
            export_svg(name, curve, points)
        # print(name + ix_id)

        # err = curve_error_max(points, curve)

        '''
        if 0:
            # simple but we can get more useful output
            self.assertEqual(ix_final, ix_naive)
        else:
            # prints intersect points where differ,
            # more useful for debugging.
            ix_final_only = tuple(sorted(set(ix_final) - set(ix_naive)))
            ix_naive_only = tuple(sorted(set(ix_naive) - set(ix_final)))
            self.assertEqual((len(ix_final), ix_final_only, ix_naive_only), (len(ix_naive), (), ()))
        '''


class FreehandTest(unittest.TestCase, TestDataFile_Helper):
    """
    Test freehand curve fitting.
    """

    def test_curve_01(self):
        self.assertTestData("test_curve_freehand_01", 0.01)
