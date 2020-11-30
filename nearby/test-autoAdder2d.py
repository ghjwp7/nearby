#!/usr/bin/python3
# -*- mode: python;  coding: utf-8 -*-

'''Some tests for circumcircle routines in autoAdder2d.py'''

from pypevue import Point
from autoAdder2d import CircumCircle2, CircumCircle3, Face
from random import seed, random
import time, sys

def makeTestData(npoints, nfaces, style='xy', salt=123457, scale=10):
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

def visData(points, faces, pIn, pOut):
    import datetime
    from math import atan2, degrees, sqrt
    dt = datetime.datetime.today().strftime('%Y-%m-%d  %H:%M:%S')
    nvo = nfo = 0;   scale = 10
    filename = 'ws.scad'
    colist = ('Black','Red','Green','Yellow','Blue','Magenta','Cyan','White','Orange')
    cocount, x03 = len(colist), Point(0.3, 0, 0)
    with open(filename, 'w') as fout:
        fout.write (f'''// File {filename}, generated  {dt}
// Number of sides for round things
$fn=31;
// Diameter of dest. end of each cylinder
cylFarEnd={scale/313:0.3f};
// Diameter of sphere at vertex
ballSize={scale/239:0.3f};
// Height of Vertex label text
textSizeV={scale/29:6.3f};
// Height of Face label text
textSizeF={scale/23:6.3f};
// Height of In/Out label text
textSizeIO={scale/41:6.3f};
// File {filename} from writeSCADcode
module oneVert(trans, colo)
   translate (v=trans) color(c=colo) sphere(d=ballSize);

module oneCyl(trans, zAngle, colo, cylLen)
   translate (v=trans) rotate(a=[0,90,zAngle]) color(c=colo)
       cylinder(d1=.1,d2=cylFarEnd, h=cylLen);

module oneLabel (trans, colo, siz, txt) 
   translate (v=trans) color(c=colo)
       linear_extrude(0.06) text(size=siz, text=txt);

''')
        for i, v in enumerate(points):
            nvo += 1
            fout.write (f'''  oneVert([{str(v)}], "{colist[i%cocount]}");\n''')
            fout.write (f'''  oneLabel([{str(v+x03)}], "{colist[i%cocount]}", textSizeV, "V{i}");\n''')
            print (f'Point {i:<2} :  {str(v)}')
        
        for i, v in enumerate(pIn):
            fout.write (f'''  oneLabel([{str(v)}], "{colist[i%cocount]}", textSizeIO/4, "+");\n''')
            fout.write (f'''  oneLabel([textSizeIO/6+{str(v)}], "{colist[i%cocount]}", textSizeIO, "{i}");\n''')

        for i, v in enumerate(pOut):
            fout.write (f'''  oneLabel([{str(v)}], "{colist[i%cocount]}", textSizeIO/4, "x");\n''')
            fout.write (f'''  oneLabel([textSizeIO/6+{str(v)}], "{colist[i%cocount]}", textSizeIO, "{i}");\n''')
            
        for i, f in enumerate(faces):
            nfo += 1
            #print (f'Face {i:<2} : {f}      {colist[i%cocount]}')
            pp = points[f.p3]      # Last corner is corner before first corner
            for vnum in f.get123:
                np = points[vnum]
                dcc = np - pp      # difference of consecutive corners
                zAngle = f'{round(degrees(atan2(dcc.y, dcc.x)), 2):6.2f}'
                #cylLen = f'{i/30+sqrt(dcc.inner(dcc)):6.3f}'
                cylLen = f'{dcc.mag():6.3f}'
                fout.write (f'''  oneCyl([{str(pp)}], {zAngle}, "{colist[i%cocount]}", {cylLen});\n''')
                #print (f'     {i:<2} : {f}  @ {str(pp)} > {zAngle}  L {cylLen}')
                pp = np
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
    points, faces, pIn, pOut = makeTestData(nverts, nfaces, style=style)
    
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
    else:                       # Visualization (generates SCAD code)
        visData(points, faces, pIn, pOut)
