from math import *
from Spheral3d import *
from PolyhedronFileUtilities import *

zoidberg = readPolyhedronOBJ("zoidberg.obj")

planes = vector_of_Plane()
nplanes = 100
dtheta = 2.0*pi/nplanes
for i in xrange(nplanes):
    theta = i*dtheta
    rhat = Vector(cos(theta), sin(theta), 0.0)
    planes.append(Plane(point = 1*rhat, normal=-rhat))

clipFacetedVolumeByPlanes(zoidberg, planes)
writePolyhedronOBJ(zoidberg, "zoidberg_clipped.obj")
