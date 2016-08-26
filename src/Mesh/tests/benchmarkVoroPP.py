from Spheral3d import *
import random, time

ngens = 50000
ngencheck = ngens
nxcell = KeyTraits.maxKey1d/4
nxycell = nxcell*nxcell
ncell = nxycell*nxcell
dx = 1.0/nxcell

# Randomly generate some generators in a unit cube.
print "Creating %i generators." % ngens
generators = vector_of_Vector()
for i in random.sample(xrange(0, ncell), ngens):
    ix = i % nxcell
    iy = (i % nxycell) // nxcell
    iz = i // nxycell
    generators.append(Vector((ix + 0.5)*dx,
                             (iy + 0.5)*dx,
                             (iz + 0.5)*dx))
print "Done."

# Try out various methods of using VoroPP to break this thing up.
for nx in (10, 20, 40, 80):
    print "Trying %i..." % nx

    t0 = time.clock()
    vpp = VoroPP(generators, Vector(0,0,0), Vector(1,1,1), nx, nx, nx)
    print "   Required %g seconds to create container." % (time.clock() - t0)

    t0 = time.clock()
    a = vector_of_vector_of_Vector()
    b = vector_of_vector_of_vector_of_unsigned()
    vpp.allCells(a, b)
    print "   Required %g seconds to call allCells." % (time.clock() - t0)
