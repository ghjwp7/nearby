#!/usr/bin/env python3
'''Tests for kNN.py and related software'''
# jiw 1 Dec 2020

__author__ = "James Waldby"
__copyright__ = "James Waldby"
__license__ = "mit"

import unittest, os
from copy import deepcopy
from base_test import BaseTest
from pypevue import Point
from nearby.kNN import Cell, PNN, doAllPairs, doAMethod
from nearby.kNNmain import makeTestData
from nearby import kNNvis

class kNN_Test(BaseTest):
    '''(In following, nearby is the project dir, not nearby/src/nearby)
      -- to run just this test set:
            cd nearby;  python3 -m unittest discover tests -p kNN_test.py
              -or-
            cd nearby/tests;  ./kNN_test.py
      -- to run all tests:
            cd nearby;  python3 setup.py test
    '''
    testPath = os.path.dirname(__file__)
    visBase  = os.path.join(testPath, 'kNNtest')
    #print(f'testPath {testPath}    visBase {visBase}')
    visSetups= ((220,), {1,3,9}, True)
    
    def test_00_instantiate(self):
        print('\nPoint, kNN PNN, and kNN Cell instantiation tests')
        p = Point(1,2,3)
        q = PNN(2,3,4)
        v = Cell(7)

    def test_20_AvsB(self):  self.doStyle(0, 2)
    def test_21_AvsB(self):  self.doStyle(1, 2)
    def test_22_AvsB(self):  self.doStyle(2, 2)
    def test_23_AvsB(self):  self.doStyle(3, 2)
    def test_24_AvsB(self):  self.doStyle(4, 2)
    def test_25_AvsB(self):  self.doStyle(5, 2)

    def test_30_AvsB(self):  self.doStyle(0, 3)
    def test_31_AvsB(self):  self.doStyle(1, 3)
    def test_32_AvsB(self):  self.doStyle(2, 3)
    def test_33_AvsB(self):  self.doStyle(3, 3)
    def test_34_AvsB(self):  self.doStyle(4, 3)
    def test_35_AvsB(self):  self.doStyle(5, 3)

    def test_v32_AvsB(self):  self.doStyle(3, 2, setups=kNN_Test.visSetups)
    def test_v42_AvsB(self):  self.doStyle(4, 2, setups=kNN_Test.visSetups)
    def test_v43_AvsB(self):  self.doStyle(4, 3, setups=kNN_Test.visSetups)
    def test_v53_AvsB(self):  self.doStyle(5, 3, setups=kNN_Test.visSetups)

    def doStyle(self, dataKind, nDim, setups=((200,400), range(1,9), False)):
        '''setups is a tuple with some point counts, some k values, and a
        flag for whether to write a visualization SCAD file.'''
        countSet, nnSet, visIt = setups
        for npoints in countSet:
            print(f'\tRun kNN tests with style {dataKind} {nDim}D data, {npoints} points, k:{list(nnSet)}')
            for PNN.kNN in nnSet:   # PNN.kNN = neighbor count
                # makeTestData parameters: npoints, dstyl, nDim=2,
                #        salt=123457, scale=1, region=inUnitSquare
                verts = makeTestData(npoints, dataKind, nDim=nDim)
                ra = doAMethod(verts)
                rb = doAllPairs(verts)
                self.assertEqual(ra, rb)
                if visIt:
                    visFile = f'{self.visBase}-{nDim}-{npoints}-{dataKind}-{PNN.kNN}'
                    visForm = f'{nDim}D  n:{npoints}  s:{dataKind}  k:{PNN.kNN}'
                    kNNvis.visData(ra, visFile, dataForm=visForm)
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
