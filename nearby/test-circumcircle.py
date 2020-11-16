#!/usr/bin/env python3

'''test circumcircle code, jiw 26 Aug 2020'''

# For a triangle in R^3, if point m is the center of the triangle's
# circumcircle, we have m as follows:
#
#            |c-a|^2 [(b-a)x(c-a)]x(b-a) + |b-a|^2 (c-a)x[(b-a)x(c-a)]
#    m = a + ---------------------------------------------------------
#                            2 | (b-a)x(c-a) |^2
#
# In the above, |v| is the length (or magnitude) of vector v and x
# represents vector cross product, ie u x v is u cross v.
from sys import argv
from math import sqrt, pi, cos, sin, asin, atan2
from pypevue import Point
from random import randint as rand
#(18, n-1)
def circumcenter(a,b,c):
    '''Calc & return circumcenter point for triangle with points a,b,c'''
    cma = c.diff(a);   camm = cma.inner(cma)
    bma = b.diff(a);   bamm = bma.inner(bma)
    baxca = bma.cross(cma)
    denom = 2*baxca.inner(baxca)
    numR = cma.cross(baxca);  numR.scale(bamm/denom)
    numL = baxca.cross(bma);  numL.scale(camm/denom)
    return a.add(numR).add(numL)

    
'''
print (f'bma {bma}      bamm  {bamm}')
print (f'cma {cma}      camm  {camm}')
print (f'baxca {baxca}   denom {denom}')
print (f'numR {numR}   numL  {numR}')
tsum1 = numR.add(numL)
tsum2 = a.add(numL)
tsum3 = a.add(numR)
tsum4 = tsum3.add(numL)
print (f'tsum1  R+L {tsum1}')
print (f'tsum2  a+L {tsum2}')
print (f'tsum3  a+R {tsum3}')
print (f'tsum4 t3+L {tsum4}')
'''

a=Point(2,-3,-4);   b=Point(2,13,4);   c=Point(12,13,7)
for j in range(9):
    a, b, c = b, c, Point(rand(-5,7), rand(-5,11), rand(-13,3))
    m = circumcenter(a,b,c)
    print (f'\na {a}   b {b}   c {c}  ->  m {m} ')
    u = m.diff(a);  umm = u.mag()
    v = m.diff(b);  vmm = v.mag()
    w = m.diff(c);  wmm = w.mag()
    #print (f'Vectors from m to a, b, c:  {u}   {v}   {w} ')
    print (f'Distances from m to a, b, c:  {umm}   {vmm}   {wmm} ')
