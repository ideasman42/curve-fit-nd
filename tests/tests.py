
"""
Example of watching a single test:
  watch -n2 "USE_SVG=1 nice -n 20 python -m unittest tests.FreehandTest.test_curve_01"
"""

import os
import sys

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
            for v0, v1 in iter_pairs(s):
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
            for p in s:
                for v in p:
                    fw('<circle cx="%.4f" cy="%.4f" r="2"/>\n' %
                    (v[0] * scale, -v[1] * scale))

            # lines
            fw('<g stroke="white" stroke-opacity="0.5" stroke-width="1">\n')
            for v0, v1, v2 in s:
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v0[0] * scale, -v0[1] * scale, v1[0] * scale, -v1[1] * scale))
                fw('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" />\n' %
                (v1[0] * scale, -v1[1] * scale, v2[0] * scale, -v2[1] * scale))
            fw('</g>\n')


            fw('</g>\n')


        fw('</svg>')


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
