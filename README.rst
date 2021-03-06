.. -*- mode: rst -*-
======
nearby
======

Delaunay & k-nearest-neighbor methods

Description
===========

In its nearby/src/nearby directory, this package has Python modules
delaunay.py and kNN.py that accept lists of points (2D or 3D) and
produce either k-nearest-neighbor or Delaunay triangulation data
structures.  In its nearby/tests/ directory, it has some basic unit
tests for Delaunay and circumcircle parts of delaunay.py, and for
nearest-neighbor methods of kNN.py.

Note, kNN.py and delaunay.py both import Point (a class that
represents 3D points and operations) from pypevue, which is available
at https://github.com/ghjwp7/pypevue.git that can be installed as is
or perhaps copied to a local file.  Note, installing pypevue should
create a Python-package link to the pypevue package and a link in
.local/bin to the pypevu program.  (Alternately, use shell scripts
new-bin-lib-links and show-bin-lib-links to make or show those links.)
Similarly, installing nearby should create a Python-package link to
package nearby, which is in directory nearby/src/nearby.

To run tests, install pypevue (or otherwise import Point), cd to the
top-level nearby/ directory and say python3 setup.py test.  Or cd to
nearby/tests/ and say ./delaunay_test.py to run Delaunay triangulation
tests and ./kNN_test.py to run k-nearest-neighbor tests.  Test
routines contain numerous asserts that if satisfied indicate the code
is good to go.  Test routines also create three dozen SCAD files
containing code that OpenSCAD can interpret to provide visualizations
of test cases.

When given 3D data and setup, the Delaunay program (method
delaunay.Triangulate) uses 3D metrics but basically 2D analysis.
Thus, the triangulations it makes for 3D data may fail to meet
Delaunay criteria unless extra post-processing is done.  Note, for 2D
data, the first-nearest-neighbor graph is necessarily a subgraph of
the Delaunay triangulation for that data.  For general 3D data sets
that condition might not be necessary or typical.

Program delaunaymain.py is a standalone demo program that generates
random data sets of specified size; runs 2D or 3D **circumcircle**
tests; and can use delaunayvis.py to produce SCAD code for
visualization of solutions using OpenSCAD.  Regarding demos, see file
yume-delaunay-test.  Note, the Delaunay triangulation tests at the
moment merely check the NN graph is a subgraph of DT.  As noted above,
that condition mightn't be necessary or typical of 3D data sets.
(Anyhow, proper 3D Delaunay programs generate tetrahedrons, not just
triangles as the present program does.)

Note, the present program may suffice for 3D data sets of surfaces
like balls or polyhedrons; but see pypevue example eg-auto-freq-6-2d
for which unadjusted DT output has some struts running up and down
instead of across, the longer instead of shorter way.

Two kNN routines are given.  The slower one is brute force and takes
O(k*n*n) time (for n points and k neighbors) in current form.  It is
suitable for cross-checking output from the faster version during
software tests and for processing small point sets (eg under 100
points).  The faster routine, in current form, takes O(n*k) time and
memory when working on random-data point sets.  Note, using NN heaps
instead of NN lists in the kNN routines could reduce the O(k)
multiplier to O(log k) but that is not an immediate goal.


Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
