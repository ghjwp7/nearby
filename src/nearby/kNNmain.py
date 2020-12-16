#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This code is a test of a nearest neighbor algorithm.  It measures
# times for finding all nearest neighbors (ie the nearest neighbor for
# each point in a set of points).  In v0, test sets are 2D or 3D with
# x,y coordinates uniformly random in a unit square. In a future
# version, it might create a variety of 2D or 3D test sets, eg with
# points uniformly random in a rectangle, in a circle, on the surface
# of a hemisphere, in a sphere, in a cube, etc, plus a few non-random
# or non-uniform point distributions.

from random import seed, random
from pypevue import Point    # x,y,z point, with numerous methods
from math import sqrt, sin, cos, asin, atan2, pi, degrees
from copy import deepcopy
from nearby.kNNvis import visData
from nearby.kNN import Cell, PNN, doAMethod, doAllPairs
#----------------------------------------------------------------
def inUnitSquare(p):            # Return True if p is in region
    return 1 >= p.x >= 0 and 1 >= p.y >= 0
#----------------------------------------------------------------
def tryFax(n):  # Find product bx*by ~ n.
    sx = bx = by = tx = ty = round(n**0.5)
    bd = abs(n - tx*tx)
    while tx-ty < 11:
        while tx*ty > n:
            ty -= 1
            if bd > abs(tx*ty-n): bd, bx, by = abs(tx*ty-n), tx, ty
        tx += 1
        if bd > abs(tx*ty-n): bd, bx, by = abs(tx*ty-n), tx, ty
    return bx, by
#----------------------------------------------------------------
def makeTestData(npoints, dstyl, nDim=2, salt=123457,
                 scale=1, region=inUnitSquare):
    '''Make npoints points of specified dimension (2 or 3)
    within a specified region.'''
    seed(a=salt+npoints)
    zf = random if nDim==3 else lambda:0
    # Compute  bx, by ~ npoints  in case we need them below
    bx, by = tryFax(npoints)
    sn = int(round(npoints**0.5))

    # Note, cases with region check might use a for-loop to create ra,
    # instead of list comprehension.
    
    if dstyl==1:                # Rectangular grid of near-squares
        ra = [PNN(scale*i/(sn+random()/13), scale*j/(sn+random()/13),
                  scale*zf()) for i in range(bx) for j in range(by)]
    elif dstyl==2:              # Rectangular grid of near-triangles
        ra = [PNN(scale*(i+(j%2)/2)/(sn+random()/13), scale*((3**0.5)/2)*j/(sn+random()/13), scale*zf()) for i in range(bx) for j in range(by)]
    elif dstyl==3:              # 2x2 square of random curves
        hx, hy, dx, dy, r, ra = 0, 0, 0.23, 0.19, 1/5, []
        for j in range(npoints):
            hx = hx+r*dx*random();  hy =hy+r*dy*random()
            if abs(hx)<1 and abs(hy)<1:
                ra.append(PNN(scale*hx, scale*hy, zf()))
            else:
                dx, dy, hx, hy = -dy, dx, hx*0.7, hy*0.7
            dx, dy = dx-dy/7, dy+dx/7 # spiral left
    elif dstyl==4:                    # hemispherical shell randoms
        ra = []
        while len(ra) < npoints:
            x, y, z = scale*random(), scale*random(), scale*zf()
            r = x*x + y*y + z*z
            if 1.1 > r > 0.9:
                ra.append(PNN(x, y, abs(z)))
    elif dstyl==5:              # hemispherical shell spirals
        ra = [];  a, d, cut = 15*pi/npoints, 1/npoints, npoints//3
        for j in range(npoints):
            #if j>cut: d = .4/npoints
            x, y = cos(a*j)*j*d, sin(a*j)*j*d
            z = sqrt(1-x*x-y*y)
            if j>cut:   x += 0.26
            if j>2*cut: y += 0.27
            ra.append(PNN(x, y, z))
    else:                       # dstyl==0? - 1x1 uniform random
        ra = [PNN(scale*random(), scale*random(),
                  scale*zf()) for i in range(npoints)]
    #if not region: print(f'Created {len(ra)} {nDim}D points in style {dstyl}')
    return ra
#----------------------------------------------------------------
if __name__ == '__main__':
    from sys import argv
    import time
    arn = 0
    arn+=1; tcode = argv[arn]      if len(argv)>arn else 'bv'
    arn+=1; PCL   = argv[arn]      if len(argv)>arn else '20'
    arn+=1; nDim  = int(argv[arn]) if len(argv)>arn else 2
    arn+=1; labls = int(argv[arn]) if len(argv)>arn else 0
    arn+=1; dstyl = int(argv[arn]) if len(argv)>arn else 2
    arn+=1; kNNi  = int(argv[arn]) if len(argv)>arn else 1
    methodset = {'a':doAMethod, 'b':doAllPairs}
    # Do set of tests for each number in Point Count List
    ptime = 0
    print ('\n')
    vdeb = [None]*4999          # for debugging in a method
    for nvertices in [int(v) for v in PCL.split()]:
        PNN.kNN = kNNi
        baseName = 'wsx'
        datapoints = makeTestData(nvertices, dstyl, nDim=nDim, region=None)
        for l in tcode:
            if l in methodset:
                baseName = f'ws{l}'
                datapoints = makeTestData(nvertices, dstyl, nDim=nDim)
                # Get number of points made (may differ from nvertices)
                nverts = len(datapoints)
                baseTime = time.time()
                methodset[l](datapoints)
                ctime = time.time() - baseTime
                if ptime==0: ptime = ctime
                print (f'Test ({nDim} {labls} {dstyl} {kNNi}) for {l} with {nverts} points: {ctime:3.6f} seconds = {ctime/ptime:5.3f} x previous')
                ptime = ctime
            if l=='v':    # Visualization: Generates SCAD code in baseName file
                visData(datapoints, baseName, makeLabels=labls,
                        dataForm=f'with nDim: {nDim} Labls: {labls} Dstyl: {dstyl} kNN: {PNN.kNN}')
