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
from math import sqrt, atan2, degrees
#----------------------------------------------------------------
class Cell:
    Big = 1e99            # Big is an upper limit on numeric data
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

    def tellClosest(self, jp, verts):
        '''Return id of point in cell that's closest to vert jp'''
        p, dlo, kqlo = verts[jp], 1e99, None
        for kq in self.vList:
            q = verts[kq]
            d2 = (p-q).mag2()   # Squared distance of p and q
            if d2 < p.BSF:
                p.BSF, p.nBSF = d2, kq
            if d2 < q.BSF:
                q.BSF, q.nBSF = d2, jp
            if d2 < dlo:
                dlo, kqlo = d2, kq
        return kq        
#----------------------------------------------------------------
class PNN(Point):
    def __init__(self, x=0, y=0, z=0):
        self.x, self.y, self.z = x, y, z
        self.BSF = 1e99       # Best measure found So Far
        self.nBSF = None      # index of best neighbor so far
#----------------------------------------------------------------
def inUnitSquare(p):            # Return True if p is in region
    return 1 >= p.x >= 0 and 1 >= p.y >= 0
#----------------------------------------------------------------
def makeTestData(npoints, ndim=2, salt=123457,
                 scale=1, region=inUnitSquare):
    '''Make npoints points of specified dimension (2 or 3)
    within a specified region.'''
    seed(a=salt+npoints)
    zf = random if ndim==3 else lambda:0
    return [PNN(scale*random(), scale*random(),
                  scale*zf()) for i in range(npoints)] 
#----------------------------------------------------------------
def visData(points, baseName, makeLabels=False,
            colorFunc=lambda n1, n2: n1):
    '''Given point data with embedded nearest neighbor numbers, write
    openSCAD code for the nearest neighbor graph to file.

    points is a list of PNN objects (Point objects plus
    nearest-neighbor info.)

    baseName is a string, used as a base file name.  Output is written to
    the file named by concatenating baseName with '.scad'.

    Named parameter makeLabels is False by default.  It controls
    whether vertex-labeling code gets generated.

    Named parameter colorFunc by default is a function that returns
    its first argument.  In general, colorFunc should return an
    integer (to select a color from a color list, modulo its length)
    or a string (a color name).  colorFunc's first argument n1 is a
    vertex number.  If its second argument n2 is negative, colorFunc
    returns a vertex-label color, else returns a color for an arrow
    from vertex n1 to n2.    '''
    import datetime
    dt = datetime.datetime.today().strftime('%Y-%m-%d  %H:%M:%S')
    filename = f'{baseName}.scad';  scale = 1
    colist = ('Black','Red','Green','Yellow','Blue','Magenta','Cyan','White','Orange')
    cocount, xoff = len(colist), Point(0.01, 0, 0)
    with open(filename, 'w') as fout:
        fout.write (f'''// File {filename}, generated  {dt}
// Number of sides for round things
$fn=31;
// Main diameter of cylinder
cylMainDiam={scale/247:0.3f};
// Arrowhead max diameter
ArrowMaxDiam={scale/153:0.3f};
// Arrowhead length
ArrowLen={scale/61:0.3f};
// Diameter of sphere at vertex
ballSize={scale/139:0.3f};
// Height of Vertex label text
textSizeV={scale/29:6.3f};

module oneVert(trans, colo)
   translate (v=trans) color(c=colo) sphere(d=ballSize);

module oneArrow(trans, zAngle, colo, cylLen)
   translate (v=trans) rotate(a=[0,90,zAngle]) color(c=colo)
      union () {'{'}
         cylinder(d=cylMainDiam, h=cylLen-ArrowLen-ballSize/2);
         translate ([0, 0, cylLen-ArrowLen-ballSize/2])
            cylinder(d1=ArrowMaxDiam, d2=0, h=ArrowLen);
      {'}'}
module oneLabel (trans, colo, siz, txt) 
    translate (v=trans) color(c=colo)
        linear_extrude(0.06) text(size=siz, text=txt);

''')
        for jp, p in enumerate(points):
            if makeLabels:
                cocode = colorFunc(jp, -1)
                if type(cocode)==int:  cocode = colist[cocode%cocount]
                fout.write (f'''  oneLabel([{str(p+xoff)}], "{cocode}", textSizeV, "V{jp}");\n''')
            kq = p.nBSF
            cocode = colorFunc(jp, kq)
            if type(cocode)==int:  cocode = colist[cocode%cocount]
            fout.write (f'''  oneVert([{str(p)}], "{cocode}");\n''')
            # Make an arrow to nearest neighbor
            q = points[kq]   # q is like a nearest neighbor
            dqp = q - p      # free vector, p to q
            zAngle = f'{round(degrees(atan2(dqp.y, dqp.x)), 2):6.2f}'
            cylLen = f'{dqp.mag():6.3f}'
            fout.write (f'''  oneArrow([{str(p)}], {zAngle}, "{cocode}", {cylLen});\n''')
#----------------------------------------------------------------
def doAllPairs(verts):    # Use O(n^2) method for verification
    nv = len(verts)
    
    for jp in range(nv):
        p = verts[jp]
        for kq in range(jp+1, nv):
            q = verts[kq]
            dpq2 = (p-q).mag2()   # Squared distance of p and q
            if dpq2 < p.BSF:
                p.BSF, p.nBSF = dpq2, kq
            if dpq2 < q.BSF:
                q.BSF, q.nBSF = dpq2, jp
#----------------------------------------------------------------
def doAMethod(verts):    # A usually-faster method than brute force
    '''Puts points in locality buckets, treats each bucket, then treats
    several shells of cells around each cell.  Is faster than brute
    force method for n>29.  Time is O(n) for x-y data, vs brute force
    O(n^2).  v1 treats up to 3 shells in the x,y plane and ignores z
    shells.  Shells are treated in order of increasing distance.  A
    later version should keep going up in shells until we have found
    at least one neighbor for every point.  v1 doesn't verify that
    completion happens, but it's highly likely to happen. 

    Note, for comments about speeding up searches by distance-ordering
    of localities of points, see section 2 of "A fast all nearest
    neighbor algorithm for applications involving large point-clouds",
    by Jagan Sankaranarayanan, H. Samet, A. Varshney, 2007 (eg at
    https://www.cs.umd.edu/~varshney/papers/jagan_CG07.pdf ) .  That
    paper gives a more systematic method (usable for kNN problem)
    than the somewhat less precise, but simpler, methods here (usable
    for 1NN problem).    '''
    def makeShellList():
        # Make distance-ordered list of neighbor cells.  Start by
        # defining criterion methods for handling neighbors at
        # cardinal-and-ordinal directions CAO = [N,S,E,W,NE,SE,SW,NW].
        # Number methods num# for # in CAO return True if the
        # neighbor's row and column are in bounds.  If a cell-number
        # test fails, go on to next cell number.  Rough-distance
        # methods (thunks with embedded ruff values) return True if
        # current block's distance from proposed neighbor block is
        # below target.  Since we order the neighbors search by
        # increasing rough distance, if a rough-distance test fails,
        # all subsequent ruff tests will fail, so on failure we can
        # break out of cells search and go on to next point in verts.
        # Fine-distance methods fin# for # in CAO return True if
        # current point's distance from neighbor's nearest edge or
        # corner distance is below target. If a fine-distance test
        # fails, go on to next cell number.
        def numNE(tx, ty): return tx <= xparts and ty <= yparts
        def numSE(tx, ty): return tx <= xparts and ty >= 0
        def numSW(tx, ty): return tx >= 0      and ty >= 0
        def numNW(tx, ty): return tx >= 0      and ty <= yparts
        def numN (tx, ty): return ty <= yparts
        def numS (tx, ty): return ty >= 0
        def numE (tx, ty): return tx <= xparts
        def numW (tx, ty): return tx >= 0

        def finS (p,c,target): return (p.y-c.vHi.y)**2 < target
        def finN (p,c,target): return (p.y-c.vLo.y)**2 < target
        def finE (p,c,target): return (p.x-c.vLo.x)**2 < target
        def finW (p,c,target): return (p.x-c.vHi.x)**2 < target
        def finNE(p,c,target): return (p.x-c.vLo.x)**2 +(p.y-c.vHi.y)**2 < target
        def finSE(p,c,target): return (p.x-c.vLo.x)**2 +(p.y-c.vLo.y)**2 < target
        def finSW(p,c,target): return (p.x-c.vHi.x)**2 +(p.y-c.vLo.y)**2 < target
        def finNW(p,c,target): return (p.x-c.vHi.x)**2 +(p.y-c.vHi.y)**2 < target
        # Make shells map for first quadrant of neighbors
        smt, smp = [], []       # Shell maps temporary & permanent
        for i in range(1,xparts+1):
            ii = i*i
            smt.append((ii, i, 0))
            for j in range(1,yparts+1):
                smt.append((ii+j*j, i, j))
        print(sorted(smt))
        for dist2, kx, ky in sorted(smt):
            rufx, rufy = max(0,kx-1)*xdelt, max(0,ky-1)*ydelt
            ruff = rufx*rufx + rufy*rufy
            ruffer = lambda rufTarg, celldd=ruff: celldd < ruffTarg
            # [Following 4-symmetry assumes xdelt == ydelt] 
            # Treat 4 cells (by symmetry) at current distance
            for jj in range(4): # kx, ky = -ky, kx gets group
                cOffset = kx*xstep + ky*ystep
                grid = 1+(kx>0)-(kx<0) + 3*(1+(ky<0)-(ky>0))
                nummer = (numNW, numN, numNE,
                          numW, None,  numE,
                          numSW, numS, numSE)[grid]
                finner = (finNW, finN, finNE,
                          finW, None,  finE,
                          finSW, finS, finSE)[grid]
                todo = (cOffset, nummer, ruffer, finner, kx, ky, ruff)
                smp.append(todo)
                kx, ky = -ky, kx
        #import inspect
        for cOffset, nummer, ruffer, finner, kx, ky, ruff in smp:
            print(f'{cOffset:4}  n {nummer.__name__:5}  f {finner.__name__:5}  k {kx:2} {ky:2}   {ruff:6.3f}')
        print (xparts, yparts, xstep, ystep, zstep)
    #----------------------------------------------------------------
    #------- debug for makeShellList ---------
    xparts = yparts = 4;  xdelt  = ydelt  = 0.25
    xstep = 1; ystep = 1+xparts; zstep = ystep*(1+yparts)
    makeShellList()
    #-----------------------------------------
    nv = len(verts)
    if nv < 3: return
    # Find min & max values on each axis
    xmax, ymax, zmax = -Cell.Big, -Cell.Big, -Cell.Big
    xmin, ymin, zmin = +Cell.Big, +Cell.Big, +Cell.Big
    for q in verts:
        xmin, xmax = min(xmin, q.x),  max(xmax, q.x)
        ymin, ymax = min(ymin, q.y),  max(ymax, q.y)
        zmin, zmax = min(zmin, q.z),  max(zmax, q.z)
    #print (f'{xmin:7.3f} {xmax:7.3f} {ymin:7.3f} {ymax:7.3f} {zmin:7.3f} {zmax:7.3f}')
    xspan, yspan, zspan = xmax-xmin, ymax-ymin, zmax-zmin
    cellCount = nv/Cell.popPerCell
    if zspan==0:
        xparts = yparts = int(sqrt(cellCount)); zparts, zspan = 0, 1
    else:
        xparts = yparts = zparts = int(cellCount**(1/3))
    if xparts < 4:  xparts = yparts = 4
    xmul, ymul, zmul = xparts/xspan, yparts/yspan, zparts/zspan
    
    #print (f'{xparts:7} {xmul:7.3f} {yparts:7} {ymul:7.3f} {zparts:7} {zmul:7.3f}')
    xstep = 1; ystep = 1+xparts; zstep = ystep*(1+yparts)
    xper, yper, zper = xstep*xmul, ystep*ymul, zstep*zmul
    #print (f'{xparts:7} {xstep:7} {yparts:7} {ystep:7} {zparts:7} {zstep:7}')
    # Make lists of points in various cells.
    cellCount = zstep*(zparts+1)
    cells = [None]*cellCount
    def calcCellNum(p):
        dx, dy, dz = p.x-xmin, p.y-ymin, p.z-zmin
        #return int(dx*xper) + int(dy*yper) + int(dz*zper)
        return int(dx*xper) + int(dy*ymul)*ystep + int(dz*zmul)*zstep

    for jp in range(nv):
        cellnum = calcCellNum(verts[jp])
        if cellnum >= cellCount: continue
        if not cells[cellnum]:
            cells[cellnum] = Cell(cellnum)
        cells[cellnum].addVert(jp, verts)
    for c in cells:
        continue
        if c:
            print (f'In cell {c.cell:2}:  {c.vList} / {c.vLo} / {c.vHi}')

    for c in cells:
        if not c: continue
        cv = c.vList
        lcv = len(cv)
        #if lcv < 2:
        #    print (f'  {c.cell}:{cv}', end=''); continue
        for jp in range(lcv):
            p = verts[cv[jp]]
            for kq in range(jp+1, lcv):
                q = verts[cv[kq]]
                d2 = (p-q).mag2()   # Squared distance of p and q
                if d2 < p.BSF:
                    p.BSF, p.nBSF = d2, cv[kq]
                if d2 < q.BSF:
                    q.BSF, q.nBSF = d2, cv[jp]
    #print ()

    def roro(n1, n2):
        cn1 = calcCellNum(verts[n1])
        if n2<0:
            if len(cells[cn1].vList)>1: return "blue"
            return "red"
        cn2 = calcCellNum(verts[n2])
        if cn1==cn2: return "blue"
        return "red"
    #visData(points, "wsA1", makeLabels=True, colorFunc=roro)
    
    # Might be better to have separate x/y/z shell thicknesses
    shellThik2 = min(xspan/xparts, yspan/yparts)**2
    for c in cells:
        if not c: continue      # Skip empty cells
        cellnum = c.cell
        for jp in c.vList:
            p = verts[jp]
            # Treat neighbor cells by distance ranks (with 2-fold symmetry)
            for kx, ky, shell2 in ((0,1,0), (1,1,0), (0,2,1), (1,2,1), (2,1,1), (2,2,1), (0,3,4), (1,3,4), (3,1,4), (2,3,4), (3,2,4), (3,3,4)):
                if shell2*shellThik2 > p.BSF:
                    break
                for jj in range(4): # kx, ky = -ky, kx is for symmetry
                    dx, dy, kx, ky = kx*xstep, ky*ystep, -ky, kx
                    nbrCN = cellnum+dx+dy
                    #print (f'cellnum {cellnum:4}  nbrCN {nbrCN:4}  kx {kx:2}  ky {ky:2}  shell2 {shell2}')
                    if nbrCN<0 or nbrCN>=cellCount or not cells[nbrCN]:
                        continue
                    nbr = cells[nbrCN]
                    mx = my = 0
                    if dx > 0:    mx = nbr.vLo.x - p.x
                    elif dx < 0:  mx = p.x - nbr.vHi.x
                    if dy > 0:    my = nbr.vLo.y - p.y
                    elif dy < 0:  my = p.y - nbr.vHi.y
                    if mx*mx < p.BSF and my*my < p.BSF:
                        kq = nbr.tellClosest(jp, verts)
    #visData(points, "wsA2", makeLabels=True, colorFunc=roro)
    #visData(points, "wsA3")
        
#----------------------------------------------------------------
if __name__ == '__main__':
    from sys import argv
    import time
    arn = 0
    arn+=1; tcode = argv[arn]      if len(argv)>arn else 'bv'
    arn+=1; PCL   = argv[arn]      if len(argv)>arn else '20'
    arn+=1; ndim  = int(argv[arn]) if len(argv)>arn else 2
    arn+=1; labls = int(argv[arn]) if len(argv)>arn else 0
    methodset = {'a':doAMethod, 'b':doAllPairs}
    # Do set of tests for each number in Point Count List
    ptime = 0
    for nverts in [int(v) for v in PCL.split()]:
        baseName, points = 'wsx', makeTestData(nverts, ndim=ndim)
        for l in tcode:
            if l in methodset:
                baseName = f'ws{l}'
                points = makeTestData(nverts, ndim=ndim)
                baseTime = time.time()
                methodset[l](points)
                ctime = time.time() - baseTime
                if ptime==0: ptime = ctime
                print (f'Test time for {l} with {nverts} points: {ctime:3.6f} seconds = {ctime/ptime:5.3f} x previous')
                ptime = ctime
            if l=='v':    # Visualization: Generates SCAD code in baseName file
                visData(points, baseName, makeLabels=labls)
       
