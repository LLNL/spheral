from Spheral2d import *
import mpi
import random
from siloMeshDump import *
from math import *

rangen = random.Random()

x0, y0 = 0.0, 0.0
x1, y1 = 1.0, 1.0
nx = 10*mpi.procs

nxproc = int(sqrt(mpi.procs))
assert nxproc*nxproc == mpi.procs

dxproc = (x1 - x0)/nxproc
dyproc = (y1 - y0)/nxproc
ixproc = mpi.rank % nxproc
iyproc = mpi.rank / nxproc
xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc)
xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc)

# Randomly seed the generators.  We choose from random cells in order
# to keep nodes from getting too close together.
nxcell = KeyTraits.maxKey1d/4
nycell = nxcell
assert nx < nxcell
ncell = nxcell*nycell
dxcell = (x1 - x0)/nxcell
dycell = (y1 - y0)/nycell
xynodes_all = []
occupiedCells = set()
for k in xrange(nx*nx):
    i = rangen.randint(0, ncell)
    while i in occupiedCells:
        i = rangen.randint(0, ncell)
    ix = i % nxcell
    iy = i / nxcell
    xynodes_all.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell))
    occupiedCells.add(i)
assert len(occupiedCells) == nx*nx
xynodes_all = mpi.bcast(xynodes_all)
xynodes = [v for v in xynodes_all if testPointInBox(v, xminproc, xmaxproc)]

# Build the set of generators.
gens = vector_of_Vector()
for x in xynodes:
    gens.append(x)

mesh = PolygonalMesh(gens,
                     xmin = Vector(x0, y0),
                     xmax = Vector(x1, y1))

siloMeshDump("random_polygonal_mesh_%idomains" % mpi.procs, mesh)
