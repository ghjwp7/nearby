.. -*- mode: rst -*-
======
nearby
======

Delaunay & k-nearest-neighbor methods

Description
===========

In its src/ directory, this package has Python modules delaunay.py and
kNN.py that accept lists of points (2D or 3D) and produce either
k-nearest-neighbor or Delaunay triangulation data structures.  In its
test/ directory, it has [at present] a few basic unit tests for parts
of delaunay.py and delaunaymain.py.  Tests for kNN will be added later.

When given 3D data and setup, the Delaunay program (method
delaunay.Triangulate) uses 3D metrics but basically 2D analysis.
Thus, the triangulations it makes for 3D data may be mediocre and
might be unusable unless extra post-processing is done.

Program delaunaymain.py is a standalone demo program that generates
random data sets of specified size; runs 2D or 3D **circumcircle**
tests; and can use delaunayvis.py to produce SCAD code for
visualization of solutions using OpenSCAD.  See file
yume-delaunay-test re demos.  This needs to be modified to run
Delaunay triangulation tests as well as just circumcircle tests.

Two kNN routines are given.  The slower one is brute force and takes
O(k*n*n) time (for n points and k neighbors) in current form.  It is
suitable for cross-checking output from the faster version during
software tests and for processing small point sets (eg under 100
points).  The faster routine, in its current form, takes O(n*k) time
and memory when working on random-data point sets.  For both kNN
routines, using NN heaps instead of NN lists could reduce the O(k)
time multiplier to O(log k) but that is not an immediate goal.



Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
