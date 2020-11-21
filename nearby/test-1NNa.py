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

    def setClosest(self, jp, verts):
        '''Fix id of point in cell that's closest to vert jp'''
        p, dlo = verts[jp], 1e99
        for kq in self.vList:
            q = verts[kq]
            d2 = (p-q).mag2()   # Squared distance of p and q
            if d2 < p.BSF:
                p.BSF, p.nBSF = d2, kq
            if d2 < q.BSF:
                q.BSF, q.nBSF = d2, jp
#----------------------------------------------------------------
class PNN(Point):
    def __init__(self, x=0, y=0, z=0, noNbr=0):
        self.x, self.y, self.z = x, y, z
        self.BSF = 1e99       # Best measure found So Far
        self.nBSF = noNbr     # index of best neighbor so far
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
def makeTestData(npoints, ndim=2, salt=123457,
                 scale=1, region=inUnitSquare):
    '''Make npoints points of specified dimension (2 or 3)
    within a specified region.'''
    seed(a=salt+npoints)
    zf = random if ndim==3 else lambda:0
    # Compute  bx, by ~ npoints  in case we need them below
    bx, by = tryFax(npoints)

    # Note, for more general cases with region check, will need for-loop
    ra = [PNN(scale*random(), scale*random(),
                  scale*zf(), i) for i in range(npoints)]
    sn = int(round(npoints**0.5))
    if 0: ra = [PNN(scale*i/(sn+random()/13), scale*j/(sn+random()/13),
                  scale*zf(), i) for i in range(bx) for j in range(by)]
    if 0: ra = [PNN(scale*(i+(j%2)/2)/(sn+random()/13), scale*((3**0.5)/2)*j/(sn+random()/13), scale*zf(), i) for i in range(bx) for j in range(by)]
    return ra
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
    force method for n > ~ 80.  Time is O(n) for x-y data, vs brute
    force method's O(n^2).  In 18 Nov version, shells in the x,y plane
    are treated in order of increasing distance, with z levels not
    tested or used.

    Note, for comments about speeding up searches by distance-ordering
    of localities of points, see section 2 of "A fast all nearest
    neighbor algorithm for applications involving large point-clouds",
    by Jagan Sankaranarayanan, H. Samet, A. Varshney, 2007 (eg at
    https://www.cs.umd.edu/~varshney/papers/jagan_CG07.pdf ) .  That
    paper gives a more systematic method (usable for kNN problem) than
    the somewhat less precise, but simpler, methods here (ok for 1NN
    problem, and might extend to kNN).

    '''

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
        # Start by defining criterion methods for handling neighbors
        # at cardinal-and-ordinal directions CAO =
        # [N,S,E,W,NE,SE,SW,NW].  Number methods num# for # in CAO
        # return True if the neighbor's row and column are in bounds.
        # If a cell-number test fails, go on to next cell number.
        # Rough-distance methods (thunks with embedded ruff values)
        # return True if current block's distance from proposed
        # neighbor block is below target.  Since we order the
        # neighbors search by increasing rough distance, if a
        # rough-distance test fails, all subsequent ruff tests will
        # fail, so on failure we can break out of cells search and go
        # on to next point in verts.  Fine-distance methods fin# for #
        # in CAO return True if current point's distance from
        # neighbor's nearest edge or corner distance is below
        # target. If a fine-distance test fails, go on to next cell
        # number.
        def numNE(tx, ty): return tx < xparts and ty < yparts
        def numSE(tx, ty): return tx < xparts and ty >= 0
        def numSW(tx, ty): return tx >= 0     and ty >= 0
        def numNW(tx, ty): return tx >= 0     and ty < yparts
        def numN (tx, ty): return ty < yparts
        def numS (tx, ty): return ty >= 0
        def numE (tx, ty): return tx < xparts
        def numW (tx, ty): return tx >= 0

        def finNW(p,c,target): return (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 < target
        def finN (p,c,target): return (p.y-c.vLo.y)**2                   < target
        def finNE(p,c,target): return (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 < target
        def finE (p,c,target): return (p.x-c.vLo.x)**2                   < target
        def finSE(p,c,target): return (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2 < target
        def finS (p,c,target): return (p.y-c.vHi.y)**2                   < target
        def finSW(p,c,target): return (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 < target
        def finW (p,c,target): return (p.x-c.vHi.x)**2                   < target

        # Make shells map for first quadrant of neighbors
        smt, smp = [], []       # Shell maps temporary & permanent
        # Create temporary map to sort into shells order
        for i in range(1,xparts):
            ii = i*i
            smt.append((ii, i, 0))
            for j in range(1,yparts):
                smt.append((ii+j*j, i, j))
        # Create permanent map with 4x entries using symmetry
        for dist2, kx, ky in sorted(smt):
            rufx, rufy = max(0,kx-1)*blokSide, max(0,ky-1)*blokSide
            ruff = rufx*rufx + rufy*rufy
            for jj in range(4): # kx, ky = -ky, kx gets group
                cOffset = kx*xstep + ky*ystep
                grid = 1+(kx>0)-(kx<0) + 3*(1+(ky<0)-(ky>0))
                nummer = (numNW, numN, numNE,
                          numW, None,  numE,
                          numSW, numS, numSE)[grid]
                finner = (finNW, finN, finNE,
                          finW, None,  finE,
                          finSW, finS, finSE)[grid]
                todo = (cOffset, nummer, finner, kx, ky, ruff)
                if abs(kx) < xparts and abs(ky) < yparts:
                    smp.append(todo)
                kx, ky = -ky, kx
        return smp

    #-----------------------------------------
    # doAMethod() wants at least 2 vertices
    nv = len(verts)
    if nv < 2: return
    
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
    cellCount = nv//Cell.popPerCell
    if zspan > 0:
        vol = xspan*yspan*zspan            # occupied volume
        blokSide = (vol/cellCount)**(1/3)  # nominal cube-side
    else:
        vol = xspan*yspan                  # occupied area
        blokSide = (vol/cellCount)**(1/2)  # nominal square-side
    nx, ny, nz = [max(2,int(round(1.001*l/blokSide))) for l in (xspan,yspan,zspan)]
    if zspan == 0: nz = 1
    print (f'cellCount {cellCount}   nxyz {nx*ny*nz}   blokSide {blokSide:7.4f}')
    print ('Raw counts, spans, and extents')
    for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
        print (f'n {n:4}  s {s:7.4f}  e {n*blokSide:7.4f}')
    # Blocks may be too big or too small for computed block
    # counts, so recompute blokSide, given those counts
    blokSide = 0
    for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
        blokSide = max(blokSide, 1.001*s/n)
    print ('New counts, spans, and extents')
    for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
        print (f'n {n:4}  s {s:7.4f}  e {n*blokSide:7.4f}')
    print (f'cellCount {cellCount}   nxyz {nx*ny*nz}   blokSide {blokSide:7.4f}')

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
    for jp in range(nv):
        cellnum = calcCellNum(verts[jp]) # Put each vertex into a cell
        if cellnum >= cellCount: continue
        if not cells[cellnum]:
             # Init cell when the first vertex in it occurs
            cells[cellnum] = Cell(cellnum)
        cells[cellnum].addVert(jp, verts)
    for c in cells:
        continue # comment this line to print each cell's list & limits
        if c:
            print (f'In cell {c.cell:2}:  {c.vList} / {c.vLo} / {c.vHi}')

    #-----------------------------------------
    # Within each cell, test distances for all pairs.  For random
    # data, time per cell is O(1) on average because O(popPerCell^2)
    # is O(1), hence O(n) to do all cells on average in random case.
    # But if some cells have O(n) points and cost O(n^2), we will have
    # O(n^2) overall cost.  Generally, if O(n^p) cells have O(n^(1-p))
    # points each, cost is O(n^(2-p)), which is asymptotic to O(n^2)
    # as p goes to zero.
    for c in cells:
        if not c: continue
        cv = c.vList
        lcv = len(cv)
        for jp in range(lcv):   # For each vertex in cell, test its
            p = verts[cv[jp]]   #   distance to other vertices in cell
            for kq in range(jp+1, lcv):
                q = verts[cv[kq]]
                d2 = (p-q).mag2()   # Squared distance of p and q
                if d2 < p.BSF:      # Is q an improvement for p?
                    p.BSF, p.nBSF = d2, cv[kq]
                if d2 < q.BSF:      # Is p an improvement for q?
                    q.BSF, q.nBSF = d2, cv[jp]
    #visData(points, "wsA1", makeLabels=True, colorFunc=roro)

    #-----------------------------------------
    # Process layers or shells of cells, working outwards
    shellThik2 = min(xspan/xparts, yspan/yparts)**2
    smp = makeShellList()       # create ordered list of cells for visits
    #cshell = tuple((kx, ky, max(0,kx-1)**2+max(0,ky-1)**2) for os, nu, ru, fi, kx, ky, rf in smp)
    for c in cells:
        if not c: continue      # Skip empty cells
        cellnum = c.cell
        atx, aty = cellnum % ystep, cellnum // ystep
        for jp in c.vList:
            p = verts[jp]
            # (if p distances**2 to cell edge > p.BSF no need to
            # access nbr but we don't test that in early version)
            # Treat neighbor cells by distance ranks (with 2-fold symmetry)
            #for kx, ky, shell2 in cshell:
            for cOffset, nummer, finner, kx, ky, ruff in smp:
                if not nummer(atx+kx,aty+ky):
                    continue
                if ruff > p.BSF:
                    break
                nbrCN = cellnum+cOffset
                nbr = cells[nbrCN]
                if not nbr:
                    continue
                mx = my = 0
                if kx > 0:    mx = nbr.vLo.x - p.x
                elif kx < 0:  mx = p.x - nbr.vHi.x
                if ky > 0:    my = nbr.vLo.y - p.y
                elif ky < 0:  my = p.y - nbr.vHi.y
                if mx*mx < p.BSF and my*my < p.BSF:
                    nbr.setClosest(jp, verts)
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
                nv = len(points) # nv can differ from nverts
                baseTime = time.time()
                methodset[l](points)
                ctime = time.time() - baseTime
                if ptime==0: ptime = ctime
                print (f'Test time for {l} with {nv} points: {ctime:3.6f} seconds = {ctime/ptime:5.3f} x previous')
                ptime = ctime
            if l=='v':    # Visualization: Generates SCAD code in baseName file
                visData(points, baseName, makeLabels=labls)

