from math import *
from Spheral3d import *
from PolyhedronFileUtilities import *
from timeit import default_timer as timer

zoidberg = readPolyhedronOBJ("zoidberg.obj")
zoidbergPC = PolyClipper.Polyhedron()
PolyClipper.convertToPolyhedron(zoidbergPC, zoidberg)

# planes = vector_of_Plane()
# nplanes = 10
# dtheta = 2.0*pi/nplanes
# for i in xrange(nplanes):
#     theta = i*dtheta
#     rhat = Vector(cos(theta), sin(theta), 0.0)
#     planes.append(Plane(point = 1*rhat, normal=-rhat))


planes = vector_of_PolyClipperPlane()
planes.append(PolyClipper.PolyClipperPlane3d(Vector(0, 0, 1),
                                             Vector(0, 0, 1)))

t0 = timer()
PolyClipper.clipPolyhedron(zoidbergPC, planes)
t1 = timer()
print "PolyClipper.clipPolyhedron required", t1 - t0
chunk = Polyhedron()
PolyClipper.convertFromPolyhedron(chunk, zoidbergPC)
writePolyhedronOBJ(chunk, "zoidberg_clipped_PolyClipper.obj")

#EasyProfilerDump("clipZoidberg_timings")
Timer.TimerSummary()

# t0 = timer()
# chunk = clipFacetedVolume(zoidberg, planes)
# t1 = timer()
# print "clipFacetedVolume required", t1 - t0
# writePolyhedronOBJ(chunk, "zoidberg_clipped_r3d.obj")

