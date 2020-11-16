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
from math import exp, log, sqrt, atan2, degrees
#----------------------------------------------------------------
class NullCell:
    @property
    def cell(self): return None
    @property
    def vList(self): return []
#----------------------------------------------------------------
class Cell:
    Big = 1e99
    def __init__(self, num):
        self.cell  = num
        self.vList = []
        # Init range of coords of items in cell
        self.spanHi = Point(-Cell.Big,-Cell.Big,-Cell.Big)
        self.spanLo = Point(+Cell.Big,+Cell.Big,+Cell.Big)

    def addVert(self, jp, verts):
        self.vList.append(jp)
        # See if vertex jp increases range of coords in cell
        h, l, p = self.spanHi, self.spanLo, verts[jp]
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
    nv, Big = len(verts), 1e99
    if nv < 3: return
    # Find min & max values on each axis, and init. best-so-far data
    kq = nv-1;  q = verts[kq]
    xmin = xmax = q.x;  ymin = ymax = q.y;  zmin = zmax = q.z
    for jp in range(nv):
        p = verts[jp]
        d2 = (p-q).mag2()   # Squared distance of p and q
        #p.BSF, p.nBSF = d2, kq
        #if d2 < q.BSF:
        #    q.BSF, q.nBSF = d2, jp
        q, kq = p, jp
        xmin, xmax = min(xmin, q.x),  max(xmax, q.x)
        ymin, ymax = min(ymin, q.y),  max(ymax, q.y)
        zmin, zmax = min(zmin, q.z),  max(zmax, q.z)
    #print (f'{xmin:7.3f} {xmax:7.3f} {ymin:7.3f} {ymax:7.3f} {zmin:7.3f} {zmax:7.3f}')
    xspan, yspan, zspan = xmax-xmin, ymax-ymin, zmax-zmin
    popPerCell = 3              # desired count per cell of partition
    cellCount = nv/popPerCell
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
            print (f'In cell {c.cell:2}:  {c.vList} / {c.spanLo} / {c.spanHi}')

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

    shellThik2 = min(xspan/xparts, yspan/yparts)**2
    for c in cells:
        if not c: continue
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
                    if dx > 0:    mx = nbr.spanLo.x - p.x
                    elif dx < 0:  mx = p.x - nbr.spanHi.x
                    if dy > 0:    my = nbr.spanLo.y - p.y
                    elif dy < 0:  my = p.y - nbr.spanHi.y
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
       
