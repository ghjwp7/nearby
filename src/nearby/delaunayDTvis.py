#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# J Waldby - 3 Dec 2020

# Method visData in this module generates SCAD code for OpenSCAD
# input, to allow visualization of nearest-neighbor results. 

from pypevue import Point    # x,y,z point, with numerous methods
from math import asin, atan2, degrees
from nearby.kNN import PNN
#----------------------------------------------------------------
def visData(verts, faces, baseName, makeLabels=False, dateAttrib=True,
            dataForm='', colorFunc=lambda n1, n2: n1):
    '''Given point and face data from delaunay.triangle, write
    triangulation visualization openSCAD code to file.

    verts is a list of Vert objects (Point objects plus id number)

    faces is a list of Face objects (ordered triples of vert numbers)

    baseName is a string, used as a base file name.  Output is written to
    the file named by concatenating baseName with '.scad'.

    Named parameter makeLabels is False by default.  It controls
    whether vertex-labeling code gets generated.

    Named parameter dateAttrib is True by default.  It controls
    whether a date (of form like 2020-12-02 20:31:30) is part of line
    1 of output (an attributes-documentation line).

    Named parameter dataForm is a string included verbatim as part of
    line 1 of output.  Typically, use f'with nDim: {nDim} Labls:
    {labls} Dstyl: {dstyl} kNN: {PNN.kNN}' to generate dataForm, to
    document parameters that were used in a makeTestData() call.
    
    Named parameter colorFunc by default is a function that returns
    its first argument.  In general, colorFunc should return an
    integer (to select a color from a color list, modulo its length)
    or a string (a color name).  colorFunc's first argument n1 is a
    vertex number.  If its second argument n2 is negative, colorFunc
    returns a vertex-label color, else returns a color for an arrow
    from vertex n1 to n2.    '''
    import datetime
    dt = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S') if dateAttrib else ''
    filename = f'{baseName}.scad'
    scale = 1
    scale = 9/len(verts)**0.4   # Ad hoc formula; smaller scale for larger n
    colist = ('Black','Red','Green','Yellow','Blue','Magenta','Cyan','White','Orange')
    cocount, xoff = len(colist), Point(0.01, 0, 0)
    with open(filename, 'w') as fout:
        fout.write (f'''// File {filename}, generated {dt} by kNN {dataForm}
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

module oneArrow(trans, yAngle, zAngle, colo, cylLen, head, tail)
   translate (v=trans) rotate(a=[0,yAngle,zAngle]) color(c=colo)
      union () {'{'}
         cylinder(d=cylMainDiam, h=cylLen-ArrowLen-ballSize/2);
         translate ([0, 0, cylLen-ArrowLen-ballSize/2])
            cylinder(d1=ArrowMaxDiam, d2=0, h=ArrowLen);
      {'}'}
module oneLabel (trans, colo, siz, txt)
    translate (v=trans) color(c=colo)
        linear_extrude(0.006) text(size=siz, text=txt);

''')
        for p in verts:
            jp = p.num
            if makeLabels:
                cocode = colorFunc(jp, -1)
                if type(cocode)==int:  cocode = colist[cocode%cocount]
                fout.write (f'''  oneLabel([{str(p+xoff)}], "{cocode}", textSizeV, "V{jp}");\n''')
            cocode = colorFunc(jp, None)
            if type(cocode)==int:  cocode = colist[cocode%cocount]
            fout.write (f'''  oneVert([{p}], "{cocode}");\n''')

        for nf, f in enumerate(faces): # Make arrows to verts of faces
            cornerNums = f.get123      # Get triple of vert indices
            p = verts[cornerNums[-1]]
            cocode = colorFunc(nf, -2)
            if type(cocode)==int:  cocode = colist[cocode%cocount]
            for kq in cornerNums:
                q = verts[kq]
                dqp = q - p     # free vector, p to q
                L = dqp.mag()
                print (f'nf={nf} cornerNums={cornerNums} p={p}\t q={q}\t q#{q.num} kq={kq} L={L:0.3}')
                yAngle = f'{round(90-degrees(asin(min(1, max(-1, dqp.z/L)))), 2)}'
                zAngle = f'{round(degrees(atan2(dqp.y, dqp.x)), 2)}'
                fout.write (f'''  oneArrow([{str(p)}], {yAngle}, {zAngle}, "{cocode}", {L:6.3f}, {jp}, {kq});\n''')
                p = q
#----------------------------------------------------------------
