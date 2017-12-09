from Spheral3d import *
from PolyhedronFileUtilities import *

zoidberg = readPolyhedronOBJ("zoidberg.obj")

plane1 = Plane(point=Vector(0,0,2), normal=Vector(0,0,1))
plane2 = Plane(point=Vector(0,0,2), normal=Vector(0,0,-1))

planes = vector_of_Plane(1, plane1)
print planes
chunk1 = Polyhedron(zoidberg)
clipFacetedVolumeByPlanes(chunk1, planes)

planes = vector_of_Plane(1, plane2)
print planes
chunk2 = Polyhedron(zoidberg)
clipFacetedVolumeByPlanes(chunk2, planes)

writePolyhedronOBJ(chunk1, "zoidberg_chunkONE.obj")
writePolyhedronOBJ(chunk2, "zoidberg_chunkTWO.obj")
