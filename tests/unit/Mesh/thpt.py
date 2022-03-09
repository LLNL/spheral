from Spheral3d import *
from siloMeshDump import *
import random

eos = GammaLawGasMKS(5.0/3.0, 1.0)

def createNodes(gens):
    nodes = makeFluidNodeList("some_nodes", eos,
                              numInternal = len(gens))
    pos = nodes.positions()
    H = nodes.Hfield()
    mass = nodes.mass()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    for i in xrange(len(gens)):
        xi = gens[i]
        pos[i] = xi
        H[i] = SymTensor.one
        mass[i] = 1.0
        rho[i] = 1.0 + xi.magnitude2()
        vel[i] = xi
    return nodes

generators = vector_of_Vector()
nx = 4
dx = 1.0
print "Creating generators for regular mesh."
for iz in xrange(nx):
    for iy in xrange(nx):
        for ix in xrange(nx):
            generators.append(Vector((ix + 0.5)*dx,
                                     (iy + 0.5)*dx,
                                     (iz + 0.5)*dx))
nodes = createNodes(generators)
print "Generating regular mesh."
mesh = PolyhedralMesh(generators, Vector(0,0,0), Vector(nx,nx,nx))
print "Writing..."
siloMeshDump("testPolyhedralHexes", mesh, nodeLists = [nodes])
print "Done."
del nodes

generators = vector_of_Vector()
n = nx**3
print "Creating generators for random mesh."
rangen = random.Random()
nxcell = KeyTraits.maxKey1d/4
nxycell = nxcell**2
assert nx < nxcell
ncell = nxcell**3
dxcell = 1.0/nxcell
dycell = 1.0/nxcell
dzcell = 1.0/nxcell
occupiedCells = set()
for i in xrange(n):
    i = rangen.randint(0, ncell)
    while i in occupiedCells:
        i = rangen.randint(0, ncell)
    ix = i % nxcell
    iy = (i % nxycell) / nxcell
    iz = i / nxycell
    generators.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell, (iz + 0.5)*dzcell))
    occupiedCells.add(i)
assert len(occupiedCells) == n
nodes = createNodes(generators)
print "Generating random mesh."
mesh = PolyhedralMesh(generators, Vector(0,0,0), Vector(1,1,1))
print "Writing..."
siloMeshDump("testPolyhedralRandom", mesh, nodeLists = [nodes])
print "Done."

