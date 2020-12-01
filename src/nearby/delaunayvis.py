# -*- mode: python;  coding: utf-8 -*-


from pypevue import Point
from delaunay import CircumCircle2, CircumCircle3, Face

def visData(points, faces, pIn, pOut):
    '''Visualization for output from delaunay.py: write openscad scad code
corresponding to vertex and triangle data in `points` and `faces`.
The data need not be a Delaunay triangulation; for example, it could
instead be CircumCircleX() test data.  Params: `points`, vertices
data; `faces`, triangle data; `pIn`, list of Points to label with +;
`pOut`, list of Points to label with x.  Note, tests generate pIn &
pOut inside or outside each face, not as vertices of a face.    '''
    import datetime
    from math import atan2, degrees, sqrt
    dt = datetime.datetime.today().strftime('%Y-%m-%d  %H:%M:%S')
    nvo = nfo = 0;   scale = 10
    filename = 'ws.scad'        # output-file name
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
                np = points[vnum]
                dcc = np - pp      # difference of consecutive corners
                zAngle = f'{round(degrees(atan2(dcc.y, dcc.x)), 2):6.2f}'
                cylLen = f'{dcc.mag():6.3f}'
                fout.write (f'''  oneCyl([{str(pp)}], {zAngle}, "{colist[i%cocount]}", {cylLen});\n''')
                pp = np
#--------------------------------------------------------------
