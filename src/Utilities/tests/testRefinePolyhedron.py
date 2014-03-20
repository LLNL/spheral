import time
from Spheral3d import *

def PLC_from_Polyhedron(poly):
    coords = vector_of_double()
    for v in poly.vertices():
        coords.append(v.x)
        coords.append(v.y)
        coords.append(v.z)
    plc = polytope.PLC3d()
    for f in poly.facets():
        indices = f.ipoints
        plc.facets.push_back(vector_of_int())
        for i in indices:
            plc.facets[-1].append(i)
    assert len(plc.facets) == len(poly.facets())
    return coords, plc

# Create a simple unit box.
verts0 = vector_of_Vector()
verts0.append(Vector(0,0,0))
verts0.append(Vector(1,0,0))
verts0.append(Vector(1,1,0))
verts0.append(Vector(0,1,0))
verts0.append(Vector(0,0,1))
verts0.append(Vector(1,0,1))
verts0.append(Vector(1,1,1))
verts0.append(Vector(0,1,1))
poly0 = Polyhedron(verts0)

# Now refine it.
t0 = time.clock()
poly1 = refinePolyhedron(poly0, 5)
t1 = time.clock()
print "Requred %g seconds to refine polyhedron." % (t1 - t0)
coords1, plc1 = PLC_from_Polyhedron(poly1)
polytope.writePLCtoOFF(plc1, coords1, "refine_cube.off")


