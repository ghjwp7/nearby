#!/usr/bin/env python3
# -*- coding: utf-8 -*- -- jiw 15 Oct 2020

# In this module:

# Triangulate(pxy), a Delaunay triangulation routine.  It modifies pxy in place 

# Face: a class to contain points of a triangular face

# Vert:  a class to contain 1 point and its index number

# CircumCircle2(point, threepoints, canon, cache),
# CircumCircle3(point, threepoints, canon, cache):
#   CC2 and CC3 are functions to test if `point` is in the
#   circumcircle of points indexed by `threepoints`.  `canon` is a
#   canonical signature for the three points, and `cache` is a dict
#   with circumcircle data for previously-encountered point triples.

from math import sqrt
from pypevue import Point
#==============================================================
class Face:                     # Class for triangular faces
    def __init__(self, p1,p2,p3):
        self.p1, self.p2, self.p3 = p1, p2, p3
        self.complete = False
    @property    # t.get123 is the point numbers of triangle t
    def get123(self): return self.p1, self.p2, self.p3
    def __str__(t):   return f'({t.p1:3}, {t.p2:3}, {t.p3:3})'
    def origNum(self, A):
        p,q,r = self.get123
        u = A[p].num if type(A[p])==Vert else -1
        v = A[q].num if type(A[q])==Vert else -1
        w = A[r].num if type(A[r])==Vert else -1
        return f'({p:3}, {q:3}, {r:3}) : [{u:3}, {v:3}, {w:3}]'
    @property    # t.canon is a canonical signature for point#s of t
    def canon(self):
        s = self.get123
        return (min(s), max(s), sum(s))
#==============================================================
class Vert(Point):
    def __init__(self, p, num):
        self.x, self.y, self.z, self.num = p.x, p.y, p.z, num
#==============================================================
def Triangulate(pxy):
    '''Accepts a list of vertices in array pxy.  Returns a list of
    triangular faces with vertices going clockwise. '''
    #  Note, Triangulate mods pxy by sorting it and by adding 3
    # `super-triangle` points at end
    
    #  Returns a list of triangular Face elements, with vertices going
    #  `clockwise` which may be ambiguous for folded 3D structures.  Do
    #  your own orientation if it matters.
    
    # The code in this half-3D routine was adapted from
    # triangulate.py, a 2D 22 Oct 2020 python version by James Waldby
    # of code in 2D Delauney triangulation programs by: Paul Bourke,
    # c, 1989; fjenett, java, 2005; Gregory Seidman, ruby, 2006.  See
    # links at http://paulbourke.net/papers/triangulate/

    pxy.sort(key=lambda p: p.x) # Sort points on ascending x coordinate
    # Allocate memory for the completeness list, flag for each triangle
    nv = len(pxy)
    #  Find maximum and minimum vertex bounds, for calculation of
    #  bounding triangle (supertriangle)
    xmin, xmax = min(v.x for v in pxy), max(v.x for v in pxy)
    ymin, ymax = min(v.y for v in pxy), max(v.y for v in pxy)
    dx,  dy = xmax - xmin, ymax - ymin
    dmax = dx if dx > dy else dy
    xmid, ymid = (xmax + xmin) / 2.0, (ymax + ymin) / 2.0

    #print (f'\n xmin:{xmin:0.2f}   xmax:{xmax:0.2f}   ymin:{ymin:0.2f}   ymax:{ymax:0.2f}  \n dx:{dx:0.2f}   dy:{dy:0.2f}   dmax:{dmax:0.2f}   xmid:{xmid:0.2f}   ymid:{ymid:0.2f}\n')
    
    #  Set up the supertriangle, a triangle to encompass all the sample
    #  points.  Its coordinates are added at the end of the vertex list
    #  and as the first triangle in the triangle list.
    pxy.append(Point(xmid - 20 * dmax, ymid - dmax))
    pxy.append(Point(xmid,        ymid + 20 * dmax))
    pxy.append(Point(xmid + 20 * dmax, ymid - dmax))
    print (f'Super-triangle: ({pxy[-3]}), ({pxy[-2]}), ({pxy[-1]})\n')
    tris = [Face(nv, nv+1, nv+2)] # Start tris list with super-triangle
    cache = {}
    #  Add points one by one into tris, the present state of mesh
    for i in range(nv):
        p = pxy[i]              # Get next point to work on Set up
        #  edge list.  If point p is inside circumcircle of some
        #  triangle T, then edges of T get put into edge list, and T
        #  gets removed from tris.
        j, edges = 0, []  # edges are in class IEDGE
        while j < len(tris):
            if tris[j].complete:
                j += 1; continue # Skip tests if triangle j already done
            threep = [pxy[p] for p in tris[j].get123]
            # Get is-inside flag, center point c, r*r, d*d
            inside, c, rr, dd = CircumCircle3(p, threep, tris[j].canon, cache)
            #print (f'inside:{inside}   c:{c}   rr:{rr}   dd:{dd}')
            if c.x < p.x and dd > rr:
                tris[j].complete = True # Due to x-sorting of pxy, j is done
            if inside:             # Add 3 edges into edge list
                p1, p2, p3 = tris[j].get123
                edges.append((p1, p2))
                edges.append((p2, p3))
                edges.append((p3, p1))
                tri = tris.pop()
                if len(tris) <= j: break
                tris[j] = tri
                continue        # don't increase j
            j += 1

        # Tag [ie mark out] multiple edges [edges shared by different
        # triangles].  Note: if all triangles are specified
        # anticlockwise then all interior edges are opposite pointing
        # in direction, as tested for first.
        nedge = len(edges)
        for j in range(nedge-1):
            for k in range(j+1, nedge):
                # 1-2, 2-1 case ...
                if edges[j][0]==edges[k][1] and edges[j][1]==edges[k][0]:
                    edges[j] = edges[k] = (-1,-1) # Remove jth & kth edges

                #  1-1, 2-2 case ... shouldn't need it, see above
                if edges[j][0]==edges[k][0] and edges[j][1]==edges[k][1]:
                    edges[j] = edges[k] = (-1,-1) # Remove jth & kth edges

        # Form new clockwise triangles for the current point, skipping
        # tagged edges.
        for j in range(nedge):
            if edges[j][0] < 0 or edges[j][1] < 0:
                continue
            tris.append(Face(edges[j][0], edges[j][1], i))

    # Remove triangles with supertriangle vertices (triangles which
    # have a vertex numbered > nv)
    o = 0
    for t in tris:
        if t.p1 < nv and t.p2 < nv and t.p3 < nv:
            tris[o] = t
            o += 1
    tris = tris[:o]
    return tris, cache
#==============================================================
def CircumCircle2(p, threep, canon, cache):
    # Returns a four-tuple, (inside, xc, yc, rr) where `inside` is
    # True if point p is inside the circumcircle that points p1,
    # p2, p3 define.  `xc, yc` is circumcircle center.  `rr` is
    # squared radius r of circumcircle.  Points on edge of circle
    # register as inside it.  xc, yc, and r are obtained from cache if
    # possible.
    EPSILON = 1e-8              # (not specified in C program??)
    xp, yp = p.x, p.y
    if canon in cache:
        cctr, rr = cache[canon]
        dd = (xp-cctr.x)**2+(yp-cctr.y)**2
        return dd-rr < EPSILON, cctr, rr, dd

    # Check for coincident points
    p1, p2, p3 = threep
    if abs(p1.y-p2.y) < EPSILON and abs(p2.y-p3.y) < EPSILON:
        return False,(0,0),0,0

    if abs(p2.y-p1.y) < EPSILON:
        m2 = - (p3.x-p2.x) / (p3.y-p2.y)
        mx2 = (p2.x + p3.x)/2
        my2 = (p2.y + p3.y)/2
        xc = (p2.x + p1.x)/2
        yc = m2 * (xc - mx2) + my2
    elif abs(p3.y-p2.y) < EPSILON:
        m1 = - (p2.x-p1.x) / (p2.y-p1.y)
        mx1 = (p1.x + p2.x)/2
        my1 = (p1.y + p2.y)/2
        xc = (p3.x + p2.x)/2
        yc = m1 * (xc - mx1) + my1
    else:
        m1 = - (p2.x-p1.x) / (p2.y-p1.y)
        m2 = - (p3.x-p2.x) / (p3.y-p2.y)
        mx1 = (p1.x + p2.x)/2
        mx2 = (p2.x + p3.x)/2
        my1 = (p1.y + p2.y)/2
        my2 = (p2.y + p3.y)/2
        xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2)
        yc = m1 * (xc - mx1) + my1

    rsqr  = (p2.x-xc)**2 + (p2.y-yc)**2
    drsqr = (xp-xc)**2   + (yp-yc)**2
    # Cache the center and rsqr
    cctr = Point(xc,yc,0)
    cache[canon] = (cctr, rsqr)
    # Return true if point is in or on circumcircle
    return drsqr-rsqr < EPSILON, cctr, rsqr, drsqr
#==============================================================
def CircumCircle3(p, threep, canon, cache):
    '''Calc circumcenter point for triangle with points threep; return
    (inside,c,rr,dd) where inside=(|p-c|<r); c=circumcenter; rr=r*r;
    r=circumcircle radius; dd=|p-c|^2.'''
    EPSILON = 1e-8              # (not specified in C program??)
    if canon in cache:
        cctr, rr = cache[canon]
        dd = (p-cctr).mag2()
        return dd-rr < EPSILON, cctr, rr, dd
    a, b, c = threep
    cma = c-a;   camm = cma*cma # For vectors, u*v is inner(u,v)
    bma = b-a;   bamm = bma*bma
    baxca = bma & cma           # For vectors, u&v is cross(u,v)
    denom = 2*(baxca.mag2())    # mag2 is magnitude squared
    numR = (bamm/denom)*(cma & baxca)
    numL = (camm/denom)*(baxca & bma)
    cctr = a + numR + numL
    rr = (b-cctr).mag2()
    dd = (p-cctr).mag2()

    rra = (a-cctr).mag2()
    rrc = (c-cctr).mag2()
    eps = rr/1e12
    if not (abs(rra-rr)<eps and abs(rrc-rr)<eps):
        print (f'cc error?  rra:{rra:0.8f}  rrb:{rr:0.8f}  rrc:{rrc:0.8f}')
    # Cache the center and rr
    cache[canon] = (cctr, rr)
    # Return true if point is in or on circumcircle
    return dd-rr < EPSILON, cctr, rr, dd
#==============================================================
