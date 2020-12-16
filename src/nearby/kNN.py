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
#----------------------------------------------------------------
class Cell:
    Big = 1e99            # Big is an upper limit on numeric data
    eps = 1e-9            # Small number for zero-testing
    popPerCell = 3        # Desired vertex count per cell of partition
    def __init__(self, num):
        self.cell  = num
        self.vList = []
        # Init range of coords of items in cell
        self.vHi = Point(-Cell.Big,-Cell.Big,-Cell.Big)
        self.vLo = Point(+Cell.Big,+Cell.Big,+Cell.Big)

    def addVert(self, jp, verts):
        self.vList.append(jp)
        # See if vertex jp increases range of coords in cell
        h, l, p = self.vHi, self.vLo, verts[jp]
        h.x, h.y, h.z = max(h.x,p.x), max(h.y,p.y), max(h.z,p.z)
        l.x, l.y, l.z = min(l.x,p.x), min(l.y,p.y), min(l.z,p.z)

    def setClosest(self, jp, verts):
        '''Fix id of point in cell that's closest to vert jp'''
        p = verts[jp]
        for kq in self.vList:
            q = verts[kq]
            d2 = (p-q).mag2()   # Squared distance of p and q
            #print (f'jp {jp}  kq {kq}  p.nBSF {p.nBSF} p.BSF {p.BSF[0]:8.5g} {p.BSF[1]:8.5g}  d2 {d2:8.5g}')
            if d2 < p.BSF[0]:  p.addNN(d2, kq)
            if d2 < q.BSF[0]:  q.addNN(d2, jp)

#----------------------------------------------------------------
class PNN(Point):
    kNN = 4               # Number of nearest neighbors we will find
    def __init__(self, x=0, y=0, z=0):
        self.x, self.y, self.z = x, y, z
        self.BSF  = [Cell.Big]*PNN.kNN  # Best measure found So Far
        self.nBSF = [-1]*PNN.kNN     # index of best neighbor so far

    def addNN(self, dist2, kq):
        # Add kq to the Best-So-Far list.  This code takes O(kNN) time
        # worst case but often takes O(1).  A heap would take
        # O(ln(kNN)) worst case.  For small kNN, O(kNN) is fast
        # enough.  Note, for random data the number of updates is
        # small because BSF[0] rapidly grows small.  Note, addNN
        # assumes dist2 is an improvement over BSF[0] -- it does not
        # check that assumption because all of the addNN() calls in
        # current code satisfy it and that way we have less overhead
        # in the kNN=1 case.  Problem: shells process app. does some
        # points twice which needs an `in` test (which can be slow) to
        # prevent duplicates on list of k closest.  Need to do that
        # test before moving any elements in array, probably.  how to
        # fix better??
        if kq in self.nBSF: return
        bat = 0
        for j in range(1, PNN.kNN):
            if self.BSF[j] > dist2:
                # If BSF[j] is superseded, move it down one level
                self.BSF[j-1], self.nBSF[j-1] = self.BSF[j], self.nBSF[j]
                bat = j
            else:
                break
        self.BSF[bat], self.nBSF[bat] = dist2, kq
#----------------------------------------------------------------
def doAllPairs(verts):    # Use O(n^2) method for verification
    nv = len(verts)

    for jp in range(nv):
        p = verts[jp]
        for kq in range(jp+1, nv):
            q = verts[kq]
            d2 = (p-q).mag2()   # Squared distance of p and q
            if d2 < p.BSF[0]:  p.addNN(d2, kq)
            if d2 < q.BSF[0]:  q.addNN(d2, jp)
    return verts
#----------------------------------------------------------------
def doAMethod(verts):    # A usually-faster method than brute force
    '''Puts points in locality buckets, treats each bucket, then treats
    several shells of cells around each cell.  Is faster than brute
    force method for n > ~ 30.  Time is O(n) for random x-y data, vs
    brute force method's O(n^2).  In current version (20 Nov), shells
    in the x,y plane are treated in order of increasing distance, with
    z levels not tested or used.

    Note, for comments about speeding up searches by distance-ordering
    of localities of points, see section 2 of "A fast all nearest
    neighbor algorithm for applications involving large point-clouds",
    by Jagan Sankaranarayanan, H. Samet, A. Varshney, 2007 (eg at
    https://www.cs.umd.edu/~varshney/papers/jagan_CG07.pdf ).  That
    paper gives a more systematic method (usable for kNN problem) than
    the somewhat less precise, but simpler, methods here (ok for 1NN
    problem, and might extend to kNN).    '''

    #-----------------------------------------
    # define color-selector function for debugging
    def roro(n1, n2):
        cn1 = calcCellNum(verts[n1])
        if n2<0:
            if len(cells[cn1].vList)>1: return "blue"
            return "red"
        cn2 = calcCellNum(verts[n2])
        if cn1==cn2: return "blue"
        return "red"
    #-----------------------------------------
    # Create distance-ordered list of neighbor cells (shells of cells)
    def makeShellList():
        # Number methods off# test if the neighbor's row and column
        # are off the grid; if so, caller goes to next cell number.
        offN = {
            -1: lambda x,y,z: y >= yparts or z <  0,
            -0: lambda x,y,z: y >= yparts,
            +1: lambda x,y,z: y >= yparts or z >= zparts}
        offNE = {
            -1: lambda x,y,z: x >= xparts or y >= yparts or z <  0,
            -0: lambda x,y,z: x >= xparts or y >= yparts,
            +1: lambda x,y,z: x >= xparts or y >= yparts or z >= zparts}
        offE = {
            -1: lambda x,y,z: x >= xparts or z <  0,
            -0: lambda x,y,z: x >= xparts,
            +1: lambda x,y,z: x >= xparts or z >= zparts}
        offSE = {
            -1: lambda x,y,z: x >= xparts or y <  0      or z <  0,
            -0: lambda x,y,z: x >= xparts or y <  0,
            +1: lambda x,y,z: x >= xparts or y <  0      or z >= zparts}
        offS = {
            -1: lambda x,y,z: y <  0      or z <  0,
            -0: lambda x,y,z: y <  0,
            +1: lambda x,y,z: y <  0      or z >= zparts}
        offSW = {
            -1: lambda x,y,z: x <  0      or y <  0      or z <  0,
            -0: lambda x,y,z: x <  0      or y <  0,
            +1: lambda x,y,z: x <  0      or y <  0      or z >= zparts}
        offW = {
            -1: lambda x,y,z: x <  0      or z <  0,
            -0: lambda x,y,z: x <  0,
            +1: lambda x,y,z: x <  0      or z >= zparts}
        offNW = {
            -1: lambda x,y,z: x <  0      or y >= yparts or z <  0,
            -0: lambda x,y,z: x <  0      or y >= yparts,
            +1: lambda x,y,z: x <  0      or y >= yparts or z >= zparts}
        offZ = {
            -1: lambda x,y,z: z <  0,
            -0: lambda x,y,z: True,
            +1: lambda x,y,z: z >= zparts }
        
        # Fine-distance methods fin# test if current point's distance
        # from neighbor's nearest datum is below target. If a
        # fine-distance test fails, caller goes to next cell number.
        finN = {
            -1: lambda p,c: (p.y-c.vLo.y)**2  + (p.z-c.vHi.z)**2  < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2                      < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2  + (p.z-c.vLo.z)**2  < p.BSF[0]}
        finNE = {
            -1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finE = {
            -1: lambda p, c:                  (p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c:                  (p.x-c.vLo.x)**2                   < p.BSF[0],
            +1: lambda p, c:                  (p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finSE = {
            -1: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2                   < p.BSF[0],
            +1: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finS = {
            -1: lambda p, c: (p.y-c.vHi.y)**2                  +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c: (p.y-c.vHi.y)**2                                    < p.BSF[0],
            +1: lambda p, c: (p.y-c.vHi.y)**2                  +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finSW = {
            -1: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finW = {
            -1: lambda p,c:                   (p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c:                   (p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c:                   (p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finNW = {
            -1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0]}
        finZ = {
            -1: lambda p,c:                                     (p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c:                                     False,
            +1: lambda p,c:                                     (p.z-c.vLo.z)**2 < p.BSF[0]}

        # Sorting function later used for de-dup because of duplicate entries
        def fDist(t):
            x, y, z = t[1], t[2], t[3]
            ss = x*x + y*y + z*z
            uid = x/100103 + y/119311 + z/139339
            return ss + uid
        
        # Make shells map for first quadrant of neighbors
        sma, smb, smc, smd = [], [], [], []    # Shell map temporaries
        # Create temporary map to sort into shells order
        for i in range(xparts):
            ii = i*i
            for j in range(yparts):
                jj = j*j
                for k in range(zparts):
                    sma.append((ii+jj+k*k, i, j, k))

        # Create real map with 4x entries using symmetry
        for dist2, kx, ky, kzin in sorted(sma):
            if kx==ky==kzin==0: continue  # all-zeroes not allowed
            # Rough distances test if current block's distance from
            # proposed neighbor block is good enough.  Neighbors
            # search is ordered by increasing ruff.  If a ruff test
            # fails, all subsequent ruff tests would also fail, so
            # break out of current point's search on first fail.
            rx, ry, rz = [max(0,k-1)*blokSide for k in (kx, ky, kzin)]
            ruff = rx*rx + ry*ry + rz*rz
            for kz in (-kzin, +kzin):
                unit = -1 if kz<0 else (1 if kz>0 else 0)
                #print (f'kx, ky, kz: {kx:3} {ky:3} {kz:3}  unit: {unit:2}  ruff: {ruff:0.4} bS {blokSide:0.4} rx {rx:0.4}  ry {ry:0.4}  rz {rz:0.4}')
                for jj in range(4): # kx, ky = -ky, kx gets group of 4
                    cOffset = kx*xstep + ky*ystep + kz*zstep
                    grid = 1+(kx>0)-(kx<0) + 3*(1+(ky<0)-(ky>0))
                    offgrid= (offNW[unit], offN[unit], offNE[unit],
                              offW [unit], offZ[unit], offE [unit],
                              offSW[unit], offS[unit], offSE[unit]) [grid]
                    finer  = (finNW[unit], finN[unit], finNE[unit],
                              finW [unit], finZ[unit], finE [unit],
                              finSW[unit], finS[unit], finSE[unit]) [grid]
                    todo = (cOffset, kx, ky, kz, ruff, offgrid, finer)
                    #print (f'kx, ky, kz: {kx:3} {ky:3} {kz:3}  unit: {unit:2}  cOffset {cOffset}  fD: {fDist(todo):0.6f}')
                    if abs(kx) < xparts and abs(ky) < yparts and abs(kz) < zparts:
                        smb.append(todo)
                    kx, ky = -ky, kx
                    if kx==ky==0: break
                if kz==0: break
        del sma; #print ()
        if 0:
            for t in smb:
                x, y, z = t[1], t[2], t[3]
                print (f'b kxyz: {x:2} {y:2} {z:2}  fD: {fDist(t):0.6f}   os: {t[0]:3}  ruff: {t[4]:0.6f}')
        if len(smb) > 0:
            smc = sorted(smb, key=fDist); smd = [smc[0]]
            del smb; #print ()
            if 0:
                for t in smc:
                    x, y, z = t[1], t[2], t[3]
                    print (f'c kxyz: {x:2} {y:2} {z:2}  fD: {fDist(t):0.6f}   os: {t[0]:3}  ruff: {t[4]:0.6f}')
            for j in range(1,len(smc)):
                if fDist(smc[j]) > fDist(smd[-1]):
                    smd.append(smc[j])
            del smc; #print ()
        else: smd = []
        if 0:
            #print (f'\n\nAfter de-dup, j={j}  len smr={len(smr)}')
            #for cOffset, kx, ky, kz, ruff, offgrid, finer in smr:
            #    print (f'fD {kx**2 + ky**2 + kz**2 - 1/cOffset:7.4f}   kxyz: {kx:3} {ky:3} {kz:3}   ruff: {ruff:0.4}')
            print (f'After de-dup, j={j}  len smd={len(smd)}  steps {ystep} {zstep}\n')
            for t in smd:
                x, y, z = t[1], t[2], t[3]
                print (f'd kxyz: {x:2} {y:2} {z:2}  fD: {fDist(t):0.6f}   os: {t[0]:3}  ruff: {t[4]:0.6f}')
        return smd

    #-----------------------------------------
    # doAMethod() wants at least 2 vertices
    nverts = len(verts)
    if nverts < 2: return
    
    #-----------------------------------------
    # Find min & max values on each axis
    xmax, ymax, zmax = -Cell.Big, -Cell.Big, -Cell.Big
    xmin, ymin, zmin = +Cell.Big, +Cell.Big, +Cell.Big
    for q in verts:
        xmin, xmax = min(xmin, q.x),  max(xmax, q.x)
        ymin, ymax = min(ymin, q.y),  max(ymax, q.y)
        zmin, zmax = min(zmin, q.z),  max(zmax, q.z)
    #print (f'{xmin:7.3f} {xmax:7.3f} {ymin:7.3f} {ymax:7.3f} {zmin:7.3f} {zmax:7.3f}')

    #-----------------------------------------
    # Compute cell block size to fit requisite blocks into occupied space
    xspan, yspan, zspan = xmax-xmin, ymax-ymin, zmax-zmin
    xflat, yflat, zflat = [(1 if abs(s)<Cell.eps else 0) for s in (xspan,yspan,zspan)]
    xmult, ymult, zmult = [(1 if abs(s)<Cell.eps else s) for s in (xspan,yspan,zspan)]
    adim = 3 - xflat - yflat - zflat
    ovol = xmult * ymult * zmult  # nominal occupied volume or area or length
    cellCount = max(1, nverts//Cell.popPerCell) # nominal cell count
    blokSide = (ovol/cellCount)**(1/adim)   # nominal cube-side
    nx, ny, nz = [max(1,int(round(1.001*l/blokSide))) for l in (xspan,yspan,zspan)]
    # Blocks may be too big or too small for computed block
    # counts, so recompute blokSide, given those counts
    blokSide = 0
    for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
        blokSide = max(blokSide, 1.001*s/n)
    if 0:
        print ('X, Y, Z  Parts, spans, and extents:')
        for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
            print (f'n {n:4}  s {s:7.4f}  e {n*blokSide:7.4f}')
        print (f'cellCount {cellCount}   nxyz {nx*ny*nz}   blokSide {blokSide:7.4f}  ^2: {blokSide**2:7.4f}')

    coMul = 1/blokSide # coMul is scale factor for x,y,z to rank,row,level
    cellCount = nx*ny*nz
    xparts, yparts, zparts = nx, ny, nz
    xstep = 1; ystep = xparts; zstep = ystep*yparts
        
    #-----------------------------------------
    # Given a point, compute its cell number
    def calcCellNum(p):
        dx, dy, dz = p.x-xmin, p.y-ymin, p.z-zmin
        return int(dx*coMul) + int(dy*coMul)*ystep + int(dz*coMul)*zstep

    #-----------------------------------------
    # Make lists of points in cells    
    cells = [None]*cellCount
    for jp in range(nverts):
        cellnum = calcCellNum(verts[jp]) # Put each vertex into a cell
        if cellnum >= cellCount: continue
        if not cells[cellnum]:
             # Init cell when the first vertex in it occurs
            cells[cellnum] = Cell(cellnum)
        cells[cellnum].addVert(jp, verts)
    if 0:
        for c in cells: # uncomment this block to print each cell's list & limits
            if c and (1 or c.cell in [2,6,7]):
                print (f'In cell {c.cell:2}:  {c.vList} / {c.vLo} / {c.vHi}')
    #-----------------------------------------
    # Within each cell, test distances for all pairs.  For random
    # data, time per cell is O(1) on average because O(popPerCell^2)
    # is O(1), hence O(n) to do all cells on average in random case.
    # But if some cells have O(n) points and cost O(n^2), we will have
    # O(n^2) overall cost.  Generally, if O(n^p) cells have O(n^(1-p))
    # points each, cost is O(n^(2-p)), which is asymptotic to O(n^2)
    # as p goes to zero.)
    for c in cells:
        if not c: continue
        cv = c.vList
        lcv = len(cv)
        for jip in range(lcv):   # For each vertex in cell, test its
            jp = cv[jip]
            p = verts[jp]
            for kiq in range(jip+1, lcv):
                kq = cv[kiq]
                q = verts[kq]
                d2 = (p-q).mag2()   # Squared distance of p and q
                if d2 < p.BSF[0]:  p.addNN(d2, kq)
                if d2 < q.BSF[0]:  q.addNN(d2, jp)
    #---debugging--------------------------------------
    def vChek(t, jp):
        for j in [0, 37, 38, 39]:
            p = verts[j]
            if vdeb[j] and vdeb[j].nBSF == p.nBSF:     continue
            vdeb[j] = deepcopy(p)
            print (f'v{j} / {jp:2} {t}  \tBSF: {[(x,round(v,3)) for x,v in zip(p.nBSF, p.BSF)]}')
    #-----------------------------------------
    # Process layers or shells of cells, working outwards
    smr = makeShellList()       # create ordered list of cells for visits
    for c in cells:
        if not c: continue      # Skip empty cells
        cellnum = c.cell
        atx, aty, atz = cellnum % ystep, (cellnum%zstep)//ystep, cellnum//zstep
        for jp in c.vList:
            p = verts[jp]
            # (if p distances**2 to cell edge > p.BSF[0] no need to
            # access nbr but we don't test that in early version)
            # Treat neighbor cells by distance ranks
            for cOffset, kx, ky, kz, ruff, offgrid, finer in smr:
                if offgrid(atx+kx, aty+ky, atz+kz):
                    continue
                #print (f'v{jp}\t{atx} {aty} {atz} : c{cellnum} to c{cellnum+cOffset}\t os {kx} {ky} {kz} : {cOffset}  \tparts: {xparts} {yparts} {zparts}   step: {ystep} {zstep}  ruff {ruff:0.4f}')
                if ruff > p.BSF[0]:
                    break
                nbr = cells[cellnum+cOffset]
                if not nbr:  continue
                #vChek('N', jp)
                mx = my = 0
                if finer(p, nbr):
                    nbr.setClosest(jp, verts)
                #vChek('F', jp)
    return verts
#----------------------------------------------------------------
if __name__ == '__main__':
    from sys import argv
    import time
    arn = 0
    arn+=1; tcode = argv[arn]      if len(argv)>arn else 'bv'
    arn+=1; PCL   = argv[arn]      if len(argv)>arn else '20'
    arn+=1; ndim  = int(argv[arn]) if len(argv)>arn else 2
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
        datapoints = makeTestData(nvertices, dstyl, ndim=ndim, region=None)
        for l in tcode:
            if l in methodset:
                baseName = f'ws{l}'
                datapoints = makeTestData(nvertices, dstyl, ndim=ndim)
                # Get number of points made (may differ from nvertices)
                nverts = len(datapoints)
                baseTime = time.time()
                methodset[l](datapoints)
                ctime = time.time() - baseTime
                if ptime==0: ptime = ctime
                print (f'Test ({ndim} {labls} {dstyl} {kNNi}) for {l} with {nverts} points: {ctime:3.6f} seconds = {ctime/ptime:5.3f} x previous')
                ptime = ctime
            if l=='v':    # Visualization: Generates SCAD code in baseName file
                visData(datapoints, baseName, makeLabels=labls)
#---------------------30---------------------
