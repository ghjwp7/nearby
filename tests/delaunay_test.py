#!/usr/bin/env python3
'''Tests for delaunay.py and related software'''
# jiw 1 Dec 2020

__author__ = __copyright__ = "James Waldby"
__license__ = "mit"

import unittest
import os, sys, random
from base_test import BaseTest
from pypevue import Point
from nearby.delaunay import Face, Vert
from nearby.delaunay import Triangulate, CircumCircle2, CircumCircle3
from nearby import delaunayDTvis, kNN, kNNmain, kNNvis

class delaunay_Test(BaseTest):
    '''(In following, nearby is the project dir, not nearby/src/nearby)
      -- to run just this test set:
            cd nearby;  python3 -m unittest discover tests -p delaunay_test.py
              -or-
            cd nearby/tests;  ./delaunay_test.py
      -- to run all tests:
            cd nearby;  python3 setup.py test
    '''

    testPath = os.path.dirname(__file__)
    visBase  = os.path.join(testPath, 'DT')
    #print(f'testPath {testPath}    visBase {visBase}')
    DTvisSetup  = (10,30,90)
    
    def xtest_00_instantiate(self):
        print('  Point and Face instantiation tests')
        p = Point(1,2,3)
        q = Face(1,2,3)
        v = Vert(p, 7)

    def test_DT20_8(self): self.doDTtest(2, 0, setup=(8,))
    #def test_DT30_8(self): self.doDTtest(3, 0, setup=(8,))
    #def test_DT21_8(self): self.doDTtest(2, 1, setup=(8,))
    #def test_DT31_8(self): self.doDTtest(3, 1, setup=(8,))
    #def test_DT21_18(self): self.doDTtest(2, 1, setup=(18,))
    #def test_DT31_18(self): self.doDTtest(3, 1, setup=(18,))
    #def test_DT20_(self): self.doDTtest(2, 0, visIt=False)
    #def test_DT30_(self): self.doDTtest(3, 0, visIt=False)

    def doDTtest(self, nDim, dataKind, setup=DTvisSetup, visIt=True):
        print(f' Test Delaunay Triangulation  setup={setup}  vis={visIt}')
        for npoints in setup:
            visFile = f'{self.visBase}-{nDim}-{dataKind}-{npoints}-DT'
            visForm = f' {nDim}D  np:{npoints}  t:{dataKind}'
            print(f'  Test {visFile}')
            points = kNNmain.makeTestData(npoints, dataKind, nDim=nDim)
            verts = [Vert(p, jp) for jp, p in enumerate(points)]
            # Make Delaunay Triangulation of points (verts)
            verts, tris, cache = Triangulate(verts)
            edges = {}
            for f in tris:            # Make a dict of tris's edges 
                cornerNums = f.get123 # Get triple of vert indices
                p = verts[cornerNums[-1]]
                jp = p.num
                for cn in cornerNums:
                    q = verts[cn]
                    kq = q.num
                    edges[(min(jp,kq), max(jp,kq))] = True
                    jp = kq
            print (f'edges={sorted(edges)}')
            # Make nearest-neighbor info 
            kNN.PNN.kNN = 1
            ra = kNN.doAMethod(points)
            # Visualize data sets now (so they exist even if errors occur)
            if visIt:
                delaunayDTvis.visData(verts, tris, visFile, makeLabels=True, dataForm=visForm)
                visFile = f'{visFile[:-3]}-NN'
                kNNvis.visData(points, visFile, makeLabels=True, dataForm=visForm)
            # Test if all the nearest neighbors appear in the DT
            for jp, p in enumerate(ra):
                kq = p.nBSF[0]
                pair = (min(jp,kq), max(jp,kq))
                print (f'Testing {jp}-{kq} via {pair} and {edges.get(pair,False)}')
                self.assertIn(pair, edges)
        
if __name__ == '__main__':
    unittest.main()
