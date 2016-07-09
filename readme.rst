
************
Curve Fit ND
************

Cubic curve fitting library.


Usage
=====

This is a high level utility library that takes a line and fits a Bézier curve to it,
in any number of dimensions - usable for freehand drawing and tracing raster images.

Intended for use cases where you simply want to input a line defined by a series of points
and get back a Bézier curve that fits within an error margin.

By supporting multiple dimensions, this allows for 2D and 3D curve fitting,
however you may want to define other properties such as radius along the curve,
color, opacity... or other properties.
Having arbitrary number of *dimensions* allow for this.


Origin
======

The method used here can be found in graphics gems ``FitCurve.c``
(by Philip J. Schneider, 1990).

This implementation was taken from OpenToonz, with some additional improvements.


Fitting Method
==============

This uses a least square solver from the original graphics gems example.

However some improvements have been made.

- Arbitrary number of dimensions.
- Replace bound-box clamping with a distance limit from the weighted center.
- Improved *fallback* methods, when the least square solution fails.

  This includes:

  - Circle fit: which accurately fits the curve to a circle. 
  - Offset fit: which uses the offset of the curve to calculate handle length.


Source Code Layout
==================

- ``c/``: the main C library.
- ``c_python_ext/``: a Python3 wrapper for the C library.
- ``tests/``: test files for the library, written in Python, using ``c_python_ext``.

