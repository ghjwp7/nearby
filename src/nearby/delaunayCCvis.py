# -*- mode: python;  coding: utf-8 -*-
# James Waldby - November 2020
from pypevue import Point
from nearby.delaunay import CircumCircle2, CircumCircle3, Face
from math import asin, atan2, degrees
import datetime

def visCCData(points, faces, pIn, pOut, dateAttrib=True,
              baseName='ws', dataForm=''):
    '''Visualization for output from delaunay.py's CircumCircle routines:
    write openscad scad code corresponding to vertex and Face triangle
    data in `points` and `faces`.  This data is for CircumCircleX()
    tests rather than for a Delaunay triangulation.  Params: • points,
    vertices data; • faces, triangle data; • pIn, list of Points to
    label with +; • pOut, list of Points to label with x; •
    dateAttrib, controls whether a date (of form like 2020-12-02
    20:31:30) is part of line 1 of output, the
    attributes-documentation line; • baseName, a string used as a base
    file name (output is written to a file named by concatenating
    baseName with '.scad'); • dataForm, string, part of line 1 of
    output. Note, makeCCtestData generated pIn & pOut inside or
    outside each face's circumcircle, not as vertices of a face.    '''
    dt = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S') if dateAttrib else ''
    nvo = nfo = 0;   scale = 10
    filename = f'{baseName}.scad' # output-file name
    colist = ('Black','Red','Green','Yellow','Blue','Magenta','Cyan','White','Orange')
    cocount, x03 = len(colist), Point(0.3, 0, 0)
    with open(filename, 'w') as fout:
        fout.write (f'''// File {filename}, generated {dt} by {dataForm}
// Number of sides for round things
$fn=31;
// Diameter of dest. end of each cylinder
cylFarEnd={scale/43:0.3f};
// Diameter of sphere at vertex
ballSize={scale/39:0.3f};
// Height of Vertex label text
textSizeV={scale/39:6.3f};
// Height of Face label text
textSizeF={scale/33:6.3f};
// Height of In/Out label text
textSizeIO={scale/41:6.3f};
// File {filename} from writeSCADcode
module oneVert(trans, colo)
   translate (v=trans) color(c=colo) sphere(d=ballSize);

module oneCyl(trans, yAngle, zAngle, colo, cylLen)
   translate (v=trans) rotate(a=[0,yAngle,zAngle]) color(c=colo)
       cylinder(d1=.1,d2=cylFarEnd, h=cylLen);

module oneLabel (trans, colo, siz, txt) 
   translate (v=trans) color(c=colo)
       linear_extrude(0.06) text(size=siz, text=txt);

''')
        for i, v in enumerate(points):
            nvo += 1
            fout.write (f'''  oneVert([{str(v)}], "{colist[i%cocount]}");\n''')
            fout.write (f'''  oneLabel([{str(v+x03)}], "{colist[i%cocount]}", textSizeV, "V{i}");\n''')
            #print (f'Point {i:<2} :  {str(v)}')

        for mark, par in (('+',pIn), ('x',pOut)):
            for i, v in enumerate(par):
                fout.write (f'''  oneLabel([{str(v)}], "{colist[i%cocount]}", textSizeIO/4, "{mark}");\n''')
                fout.write (f'''  oneLabel([textSizeIO/6+{str(v)}], "{colist[i%cocount]}", textSizeIO, "{i}");\n''')

        # For each face, draw arrows from each vertex to next vertex
        for i, f in enumerate(faces):
            nfo += 1
            # first corner = next(last corner), ie, f.p1 = next(f.p3)
            pp = points[f.p3]
            for vnum in f.get123:
                qq = points[vnum]
                dqp = qq - pp     # free vector, p to q
                L = dqp.mag()     # distance of consecutive corners
                yAngle = f'{round(90-degrees(asin(min(1, max(-1, dqp.z/L)))), 2)}'
                zAngle = f'{round(degrees(atan2(dqp.y, dqp.x)), 2)}'
                cylLen = f'{L:6.3f}'
                fout.write (f'''  oneCyl([{str(pp)}], {yAngle}, {zAngle}, "{colist[i%cocount]}", {cylLen});\n''')
                pp = qq
#--------------------------------------------------------------
