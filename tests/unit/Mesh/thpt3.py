from Spheral3d import *
from siloMeshDump import *

generators = vector_of_Vector()
nx = 2
dx = 1.0/nx
for iz in range(nx):
    for iy in range(nx):
        for ix in range(nx):
            generators.push_back(Vector((ix + 0.5)*dx,
                                        (iy + 0.5)*dx,
                                        (iz + 0.5)*dx))

mesh = PolyhedralMesh(generators, Vector(0,0,0), Vector(1,1,1))
for inode in range(mesh.numNodes):
    node = mesh.node(inode)
    print(str(node.position()), list(node.zoneIDs))

siloMeshDump("testPolyhedralHexes", mesh)
