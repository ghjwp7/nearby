#!/usr/bin/python3
# -*- mode: python;  coding: utf-8 -*-
# J Waldby - November 2020
'''Some tests for circumcircle routines in autoAdder2d.py'''

from pypevue import Point
from nearby.delaunay import CircumCircle2, CircumCircle3, Face
from nearby import delaunayCCvis
from random import seed, random
import time, sys

def makeCCTestData(npoints, nfaces, style='xy', salt=123457, scale=10):
    '''Make npoints points of specified style, xy or xyz; make nfaces
    faces, and for each face a point in & out of it.'''
    # nfaces cannot exceed (npoints choose 3)
    npC3 = npoints*(npoints-1)*(npoints-2)//6
    if nfaces > npC3:
        print (f'nfaces = {nfaces} > {npoints} C 3 = {npC3}, exiting')
        sys.exit(1)
    seed(a=salt+npoints+nfaces)
    zf = random if style=='xyz' else lambda:0
    verts = [Point(scale*random(), scale*random(),
                  scale*zf()) for i in range(npoints)]
    faces, pIn, pOut, cache = [], [], [], {}
    p2 = p3 = int(npoints*random())
    while p2==p3:
            p3 = int(npoints*random())
    for i in range(nfaces):
        while 1:
            p1, p2 = p2, p3
            while p1==p3 or p2==p3:
                p3 = int(npoints*random())
            f = Face(p1,p2,p3)
            if f.canon not in cache:
                faces.append(f)
                cache[f.canon] = 1
                break
        # Create point at barycenter of current face.  Every point
        # within the triangle is also within the circumcircle (CC).
        A,B,C = verts[p1], verts[p2], verts[p3]
        pI = 1/3*(A+B+C)

        # Create a point outside current face's CC: Because A,B,C are
        # on the CC, and pI is inside the CC, any extrapolation beyond
        # A,B,C on a line from pI is outside the CC.
        xl = 0.3;  pO = (1+xl)*B - xl*pI
        #print (f'Face {i:<2} : {f}\n   in:{pI}\n    o:{pO}')
        pIn.append(pI);  pOut.append(pO)
    return verts, faces, pIn, pOut
    
def timeTest1(CircumCircle, verts, faces, pIn, pOut, tCache, note, bmul):
    baseTime = time.time()
    cache = {}
    nface = len(faces)
    print ()
    for tn, f, pI, pO in zip(range(nface), faces, pIn, pOut):
        tnum = tn*2 + bmul*nface
        threep = [verts[p] for p in f.get123]
        if tCache:
            canon = f.canon
        else:  canon = None; cache = {}
        zin, c, rr, dd = CircumCircle(pI, threep, canon, cache)
        zot, c, rr, dd = CircumCircle(pO, threep, canon, cache)
        if not zin:
            print (f'F{tn}: Error, point {pI} showed as outside {f} = {threep}')
        if zot:
            print (f'F{tn}: Error, point {pO} showed as inside {f} = {[str(v) for v in threep]}')
    lapsTime = time.time() - baseTime
    print (f'Test time for {note}: {lapsTime:3.6f} seconds')
#--------------------------------------------------------------
if __name__ == '__main__':
    from sys import argv
    arn = 0
    arn+=1; tcode  = argv[arn]      if len(argv)>arn else 'val'
    arn+=1; nverts = int(argv[arn]) if len(argv)>arn else 5
    arn+=1; nfaces = int(argv[arn]) if len(argv)>arn else 5
    arn+=1; tcache = int(argv[arn]) if len(argv)>arn else 1
    arn+=1; haveZ  = int(argv[arn]) if len(argv)>arn else 0
    style = 'xyz' if haveZ else 'xy'
    points, faces, pIn, pOut = makeCCTestData(nverts, nfaces, style=style)
    
    def timeTest(tCache, note):
        timeTest1(CircumCircle2, points, faces, pIn, pOut, tCache, 'CC2'+note, 0)
        timeTest1(CircumCircle3, points, faces, pIn, pOut, tCache, 'CC3'+note, 2)

    def valTest(tCache):
        pass

    if tcode.startswith('tim'): # Timing test, no cache
        timeTest(False, ' time test without cache')
    elif tcode.startswith('cac'): # Timing test, with cache
        timeTest(True, ' time test using cache')
    elif tcode.startswith('val'): # Values accuracy test
        valTest(True)

    if tcode.endswith('vis'):   # Visualization (generates SCAD code)
        delaunayCCvis.visCCData(points, faces, pIn, pOut)
