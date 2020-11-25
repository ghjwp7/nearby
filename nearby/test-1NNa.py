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
def makeTestData(npoints, dstyl, ndim=2, salt=123457,
                 scale=1, region=inUnitSquare):
    '''Make npoints points of specified dimension (2 or 3)
    within a specified region.'''
    seed(a=salt+npoints)
    zf = random if ndim==3 else lambda:0
    # Compute  bx, by ~ npoints  in case we need them below
    bx, by = tryFax(npoints)
    sn = int(round(npoints**0.5))

    # Note, more-general cases with region check need for-loop
    
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
    elif dstyl==5:              # hemispherical shell spiral
        ra = [];  a, d = 15*pi/npoints, 1/npoints
        for j in range(npoints):
            x, y = cos(a*j)*j*d, sin(a*j)*j*d
            ra.append(PNN(x, y, sqrt(1-x*x-y*y)))
    else:                       # dstyl==0? - 1x1 uniform random
        ra = [PNN(scale*random(), scale*random(),
                  scale*zf()) for i in range(npoints)]
    if not region: print(f'Created {len(ra)} {ndim}D points in style {dstyl}')
    return ra
#----------------------------------------------------------------
def visData(verts, baseName, makeLabels=False,
            colorFunc=lambda n1, n2: n1):
    '''Given point data with embedded nearest neighbor numbers, write
    openSCAD code for the nearest neighbor graph to file.

    verts is a list of PNN objects (Point objects plus
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
    filename = f'{baseName}.scad'
    scale = 1
    scale = 9/len(verts)**0.4
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
textSizeV={scale/79:6.3f};

module oneVert(trans, colo)
   translate (v=trans) color(c=colo) sphere(d=ballSize);

module oneArrow(trans, yAngle, zAngle, colo, cylLen)
   translate (v=trans) rotate(a=[0,yAngle,zAngle]) color(c=colo)
      union () {'{'}
         cylinder(d=cylMainDiam, h=cylLen-ArrowLen-ballSize/2);
         translate ([0, 0, cylLen-ArrowLen-ballSize/2])
            cylinder(d1=ArrowMaxDiam, d2=0, h=ArrowLen);
      {'}'}
module oneLabel (trans, colo, siz, txt)
    translate (v=trans) color(c=colo)
        linear_extrude(0.06) text(size=siz, text=txt);

''')
        for jp, p in enumerate(verts):
            if makeLabels:
                cocode = colorFunc(jp, -1)
                if type(cocode)==int:  cocode = colist[cocode%cocount]
                fout.write (f'''  oneLabel([{str(p+xoff)}], "{cocode}", textSizeV, "V{jp}");\n''')
            kq = p.nBSF[-1]
            cocode = colorFunc(jp, kq)
            if type(cocode)==int:  cocode = colist[cocode%cocount]
            fout.write (f'''  oneVert([{str(p)}], "{cocode}");\n''')
            # Make arrows to nearest neighbors
            for kq in p.nBSF:
                if kq < 0: continue
                q = verts[kq]   # q is among nearest neighbors
                dqp = q - p      # free vector, p to q
                L = dqp.mag()
                yAngle = f'{round(90-degrees(asin(min(1, max(-1, dqp.z/L)))), 2)}'
                zAngle = f'{round(degrees(atan2(dqp.y, dqp.x)), 2)}'
                cylLen = f'{L:6.3f}'
                fout.write (f'''  oneArrow([{str(p)}], {yAngle}, {zAngle}, "{cocode}", {cylLen});\n''')

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
        
        # Fine-distance methods fin# test if current point's distance
        # from neighbor's nearest datum is below target. If a
        # fine-distance test fails, caller goes to next cell number.
        finN = {
            -1: lambda p,c: (p.y-c.vLo.y)**2  + (p.z-c.vHi.z)**2  < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2                      < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2  + (p.z-c.vLo.z)**2  < p.BSF[0] }
        finNW = {
            -1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finNE = {
            -1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vLo.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finE = {
            -1: lambda p, c:                  (p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c:                  (p.x-c.vLo.x)**2                   < p.BSF[0],
            +1: lambda p, c:                  (p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finSE = {
            -1: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2                   < p.BSF[0],
            +1: lambda p, c: (p.y-c.vHi.y)**2 +(p.x-c.vLo.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finS = {
            -1: lambda p, c: (p.y-c.vHi.y)**2                  +(p.z-c.vHi.z)**2 < p.BSF[0],
            -0: lambda p, c: (p.y-c.vHi.y)**2                                    < p.BSF[0],
            +1: lambda p, c: (p.y-c.vHi.y)**2                  +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finSW = {
            -1: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c: (p.y-c.vHi.y)**2 +(p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] }
        finW = {
            -1: lambda p,c:                   (p.x-c.vHi.x)**2 +(p.z-c.vHi.z)**2 < p.BSF[0],
            +0: lambda p,c:                   (p.x-c.vHi.x)**2 < p.BSF[0],
            +1: lambda p,c:                   (p.x-c.vHi.x)**2 +(p.z-c.vLo.z)**2 < p.BSF[0] } 

        # Make shells map for first quadrant of neighbors
        smt, smr = [], []       # Shell maps temporary & real
        # Create temporary map to sort into shells order
        for i in range(1,xparts):
            ii = i*i
            smt.append((ii, i, 0, 0))
            for j in range(1,yparts):
                jj = j*j
                smt.append((ii+jj, i, j, 0))
                for k in range(1,zparts):
                    smt.append((ii+jj+k*k, i, j, k))
                
        # Create real map with 4x entries using symmetry
        for dist2, kx, ky, kz in sorted(smt):
            # Rough distances test if current block's distance from
            # proposed neighbor block is good enough.  Neighbors
            # search is ordered by increasing ruff.  If a ruff test
            # fails, all subsequent ruff tests would also fail, so
            # break out of current point's search on first fail.
            rx, ry, rz = [max(0,k-1)*blokSide for k in (kx, ky, kz)]
            ruff = rx*rx + ry*ry + rz*rz
            for kk in (-kz, +kz):
                zz = 1 if kk else 0
                for jj in range(4): # kx, ky = -ky, kx gets group of 4
                    cOffset = kx*xstep + ky*ystep + kz*zstep
                    grid = 1+(kx>0)-(kx<0) + 3*(1+(ky<0)-(ky>0))
                    offgrid= (offNW[zz], offN[zz], offNE[zz],
                              offW [zz],   None,   offE [zz],
                              offSW[zz], offS[zz], offSE[zz]) [grid]
                    finer  = (finNW[zz], finN[zz], finNE[zz],
                              finW [zz],   None,   finE [zz],
                              finSW[zz], finS[zz], finSE[zz]) [grid]
                    todo = (cOffset, offgrid, finer, kx, ky, kz, ruff)
                    if abs(kx) < xparts and abs(ky) < yparts:
                        smr.append(todo)
                    kx, ky = -ky, kx
                if not kz: break
        return smr

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
    #nx, ny, nz = [max(1, f) for f in (xflat, yflat, zflat)]
    # Blocks may be too big or too small for computed block
    # counts, so recompute blokSide, given those counts
    blokSide = 0
    for n,s in ((nx,xspan),(ny,yspan),(nz,zspan)):
        blokSide = max(blokSide, 1.001*s/n)
    if 0:
        print ('X, Y, Z  Parts, spans, and extents:')
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
    for jp in range(nverts):
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
    #visData(verts, "wsA1", makeLabels=True, colorFunc=roro)
    #-----------------------------------------
    # Process layers or shells of cells, working outwards
    smr = makeShellList()       # create ordered list of cells for visits
    for c in cells:
        if not c: continue      # Skip empty cells
        cellnum = c.cell
        atx, aty, atz = cellnum % ystep, (cellnum//ystep)%zstep, cellnum//zstep
        for jp in c.vList:
            p = verts[jp]
            # (if p distances**2 to cell edge > p.BSF[0] no need to
            # access nbr but we don't test that in early version)
            # Treat neighbor cells by distance ranks
            for cOffset, offgrid, finer, kx, ky, kz, ruff in smr:
                if offgrid(atx+kx, aty+ky, atz+kz):
                    continue
                if ruff > p.BSF[0]:
                    break
                nbr = cells[cellnum+cOffset]
                if not nbr:
                    continue
                mx = my = 0
                if finer(p, nbr):
                    nbr.setClosest(jp, verts)
    #visData(verts, "wsA2", makeLabels=True, colorFunc=roro)
    #visData(verts, "wsA3")

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
                print (f'Test time for {l} with {nverts} points: {ctime:3.6f} seconds = {ctime/ptime:5.3f} x previous')
                ptime = ctime
            if l=='v':    # Visualization: Generates SCAD code in baseName file
                visData(datapoints, baseName, makeLabels=labls)

