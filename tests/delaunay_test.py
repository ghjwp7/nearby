#!/usr/bin/env python3
'''Tests for delaunay.py and related software'''
# jiw 1 Dec 2020

__author__ = "James Waldby"
__copyright__ = "James Waldby"
__license__ = "mit"

import unittest
import os, sys, random
from base_test import BaseTest
from pypevue import Point
from nearby.delaunaymain import makeTestData
from nearby.delaunay import Face, Vert, CircumCircle2, CircumCircle3

class delaunay_Test(BaseTest):
    '''to run:
      - cd nearby   (The project dir, not nearby/src/nearby)
      - to run just this test:
           python3 -m unittest discover tests -p test.py
      - to run all tests:
           python3 setup.py test
    '''
    testPath = os.path.dirname(__file__)
    
    def test_00_instantiate(self):
        print('\nPoint and Face instantiation tests')
        p = Point(1,2,3)
        q = Face(1,2,3)
        v = Vert(p, 7)

    def test_03_(self):
        print('\nmakeTestData - basic functionality')
        # makeTestData(npoints, nfaces, style='xy', salt=123457, scale=10)
        npoints, nfaces = 100, 2000
        # By default, data is 2D (ie z coord is zero)
        verts, faces, pIn, pOut = makeTestData(npoints, nfaces)
        self.assertEqual(npoints, len(verts))
        self.assertEqual(nfaces,  len(faces))
        self.assertEqual(nfaces,  len(pIn))
        self.assertEqual(nfaces,  len(pOut))

    def test_04_(self):
        print('\nmakeTestData + CircumCircle in/out testing')
        # makeTestData(npoints, nfaces, style='xy', salt=123457, scale=10)
        npoints, nfaces = 100, 2000
        # By default, data is 2D (ie z coord is zero)
        verts, faces, pIn, pOut = makeTestData(npoints, nfaces)
        self.assertEqual(npoints, len(verts))
        self.assertEqual(nfaces,  len(faces))
        self.assertEqual(nfaces,  len(pIn))
        self.assertEqual(nfaces,  len(pOut))
        # Tests for not-cached circumcircle calcs:
        # CircumCircle2(point, threepoints, canon, cache),
        # CircumCircle3(point, threepoints, canon, cache):
        for pIO, ioio, tfIO in ((pIn,'in', True), (pOut,'out of',False)):
            for f, p in zip(faces, pIO):
                threep = [verts[p] for p in f.get123]
                for cc in (CircumCircle2, CircumCircle3):
                    zin, ctr, rr, dd = cc(p, threep, f.canon, {})
                    self.assertEqual(zin, tfIO, f'{cc.__name__}: point {p} should be {ioio} circumcircle of face {f} = {threep}')
    '''
    def test_04_(self): pass        
    
    def test_05_(self): pass
        
    def test_06_(self): pass        
    
    def test_07_(self): pass        
    
    def test_08_(self): pass        
    
    def test_09_(self): pass        
    
    def test_10_(self): pass        
    
    def test_11_(self): pass
    '''
    
if __name__ == '__main__':
    unittest.main()
