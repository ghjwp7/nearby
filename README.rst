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
test/ directory, it has some basic unit tests for Delaunay and
circumcircle parts of delaunay.py, and for nearest-neighbor methods of
kNN.py.

When given 3D data and setup, the Delaunay program (method
delaunay.Triangulate) uses 3D metrics but basically 2D analysis.
Thus, the triangulations it makes for 3D data may fail to meet
Delaunay criteria unless extra post-processing is done.  Note, for 2D
data, the first-nearest-neighbor graph (sat 1NNg) is necessarily a
subgraph of the Delaunay triangulation for that data.  For general 3D
data sets that condition is neither necessary nor typical.

Program delaunaymain.py is a standalone demo program that generates
random data sets of specified size; runs 2D or 3D **circumcircle**
tests; and can use delaunayvis.py to produce SCAD code for
visualization of solutions using OpenSCAD.  See file
yume-delaunay-test re demos.  Note, the Delaunay triangulation tests
at the moment merely check that 1NNg is a subgraph of the DT.  As
noted above, that condition isn't necessary or typical of 3D data
sets.  (Moreover, proper 3D Delaunay programs generate tetrahedrons,
not just triangles as the present program does.)

Note, the present program may suffice for 3D data sets of surfaces
like balls or polyhedrons; but see pypevue example eg-auto-freq-6-2d
for which unadjusted DT output has some struts running up and down
instead of across, the longer instead of shorter way.

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
