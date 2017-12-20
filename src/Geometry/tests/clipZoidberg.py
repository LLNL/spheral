from math import *
from Spheral3d import *
from PolyhedronFileUtilities import *
from timeit import default_timer as timer

zoidberg = readPolyhedronOBJ("zoidberg.obj")

planes = vector_of_Plane()
nplanes = 10
dtheta = 2.0*pi/nplanes
for i in xrange(nplanes):
    theta = i*dtheta
    rhat = Vector(cos(theta), sin(theta), 0.0)
    planes.append(Plane(point = 1*rhat, normal=-rhat))

chunk = Polyhedron(zoidberg)
t0 = timer()
clipFacetedVolumeByPlanes(chunk, planes)
t1 = timer()
print "clipFacetedVolumeByPlanes required", t1 - t0
writePolyhedronOBJ(chunk, "zoidberg_clipped_native.obj")

# t0 = timer()
# chunk = clipFacetedVolume(zoidberg, planes)
# t1 = timer()
# print "clipFacetedVolume required", t1 - t0
# writePolyhedronOBJ(chunk, "zoidberg_clipped_r3d.obj")

