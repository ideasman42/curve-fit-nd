
************
Curve Fit nD
************

Cubic curve fitting library.

.. figure:: https://cloud.githubusercontent.com/assets/1869379/17393989/1f8f21ba-5a6a-11e6-82cb-9d7cda155825.png

   Example showing a cubic spline fit to 283 points.

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


Use Cases
=========

- `Blender's Freehand Drawing
  <https://archive.blender.org/wiki/index.php/Dev:Ref/Release_Notes/2.78/Modelling/#Curve_Editing>`__.
- A Rust port of this library is used in
  `Raster Retrace
  <https://github.com/ideasman42/raster-retrace>`__.


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
- Re-fitting, an alternate, a more computationally intensive method for knot placement
  which initializes the curve as a dense curve
  (where every point is a knot), then iteratively removing knots which give the least error,
  some further adjustments are made after this to avoid local-maximums giving skewed results.
  This also has the advantage that it can be used to detect corners.


Source Code Layout
==================

- ``c/``: the main C library.
- ``c_python_ext/``: a Python3 wrapper for the C library.
- ``tests/``: test files for the library, written in Python, using ``c_python_ext``.

