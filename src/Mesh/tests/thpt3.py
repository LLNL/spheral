from Spheral3d import *

generators = vector_of_Vector()
nx = 2
dx = 1.0/nx
for iz in xrange(nx):
    for iy in xrange(nx):
        for ix in xrange(nx):
            generators.push_back(Vector((ix + 0.5)*dx,
                                        (iy + 0.5)*dx,
                                        (iz + 0.5)*dx))

mesh = PolyhedralMesh(generators, Vector(0,0,0), Vector(1,1,1))
