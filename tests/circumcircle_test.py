#!/usr/bin/env python3
'''Tests for delaunay.py and related software'''
# jiw 1 Dec 2020

__author__ = __copyright__ = "James Waldby"
__license__ = "mit"

import unittest
import os, sys, random
from base_test import BaseTest
from pypevue import Point
from nearby.delaunaymain import makeCCTestData
from nearby.delaunay import Face, Vert, CircumCircle2, CircumCircle3
from nearby import delaunayCCvis

class CC_Test(BaseTest):
    '''(In following, nearby is the project dir, not nearby/src/nearby)
      -- to run just this test set:
            cd nearby;  python3 -m unittest discover tests -p delaunay_test.py
              -or-
            cd nearby/tests;  ./delaunay_test.py
      -- to run all tests:
            cd nearby;  python3 setup.py test
    '''

    testPath = os.path.dirname(__file__)
    visBase  = os.path.join(testPath, 'CC')
    #print(f'testPath {testPath}    visBase {visBase}')
    CCvisSetups  = ((10,10), (20,40))
    CCplainSetups= ((20,200), (100,2000), (200,8000))
    
    def test_00_instantiate(self):
        print('  Point and Face instantiation tests')
        p = Point(1,2,3)
        q = Face(1,2,3)
        v = Vert(p, 7)

    def test_CC03_(self):
        print(' makeTestData - basic functionality')
        # makeCCTestData(npoints, nfaces, style='xy', salt=123457, scale=10)
        npoints, nfaces = 100, 2000
        # By default, data is 2D (ie z coord is zero)
        verts, faces, pIn, pOut = makeCCTestData(npoints, nfaces)
        self.assertEqual(npoints, len(verts))
        self.assertEqual(nfaces,  len(faces))
        self.assertEqual(nfaces,  len(pIn))
        self.assertEqual(nfaces,  len(pOut))

    def doCCtest(self, dimCode, setups=CCplainSetups, visIt=False):
        '''Run a set of tests for not-cached circumcircle calcs.  For each
        Face in a set of test data, verify that a point in the
        circumcircle of the face is accepted as in the CC, and a point
        not in the CC is accepted as outside.

        Parameters: • dimCode is 'xy' or 'xyz' to signify 2D or 3D
        data.  • setups is a tuple of tuples that each contain a point
        count and a face count.  • visIt is True to cause
        visualization output, else False.      
 
        Note, each CircumCircle routine has four parameters: point,
        threepoints, canon, cache. • point is a vertex to be tested if
        in or out of CC. • threepoints is a tuple of 3 points, the
        corners of a triangular face • canon is a canonical id number
        for the three corner indices • cache is None, or a dict in
        which are cached earlier CircumCircle calcs (CC center and
        radius)        '''
        print(f' makeTestData + CircumCircle in/out testing {dimCode} data, setups={setups}  visIt={visIt}')
        for npoints, nfaces in setups:
            verts, faces, pIn, pOut = makeCCTestData(npoints, nfaces, style=dimCode)
            self.assertEqual(npoints, len(verts))
            self.assertEqual(nfaces,  len(faces))
            self.assertEqual(nfaces,  len(pIn))
            self.assertEqual(nfaces,  len(pOut))

            # Process the inside-points and outside-points lists
            for pIO, ioio, tfIO in ((pIn,'in', True), (pOut,'out of',False)):
                # For each face f and its corresponding in or out point p,
                # see if both CC routines return a status equal to tfIO.
                for f, p in zip(faces, pIO):
                    # Get the xyz coordinates of each corner of Face f 
                    threep = [verts[q] for q in f.get123]
                    for cc in (CircumCircle2, CircumCircle3):
                        # Get in or out flag (zin); CC center &
                        # radius^2, ctr & rr; and dist^2(ctr,p)
                        zin, ctr, rr, dd = cc(p, threep, f.canon, {})
                        self.assertEqual(zin, tfIO, f'{cc.__name__}: point {p} should be {ioio} circumcircle of face {f} = {threep}')
            if visIt:
                visFile = f'{self.visBase}-{dimCode}-{npoints}-{nfaces}'
                visForm = f' d:{dimCode}  np:{npoints}  nf:{nfaces}'
                delaunayCCvis.visCCData(verts, faces, pIn, pOut, baseName=visFile, dataForm=visForm)

    def test_CC04_(self): self.doCCtest('xy')
    def test_CC05_(self): self.doCCtest('xyz')
    def test_CC06v(self): self.doCCtest('xy', setups=self.CCvisSetups, visIt=True)
    def test_CC07v(self): self.doCCtest('xyz', setups=self.CCvisSetups, visIt=True)
        
if __name__ == '__main__':
    unittest.main()
