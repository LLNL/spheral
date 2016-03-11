import mpi
from math import *
from Spheral3d import *
from FacetedSurfaceGenerator import *
from siloPointmeshDump import *
from VoronoiDistributeNodes import distributeNodes3d
from PolyhedronFileUtilities import *

verts = vector_of_Vector()
for vi in [Vector(0,0,0), Vector(1,0,0), Vector(0,1,0), Vector(1,1,0), 
           Vector(0,0,1), Vector(1,0,1), Vector(0,1,1), Vector(1,1,1)]:
    verts.append(vi)
surface = Polyhedron(verts)

flags = [1, 1, 1, 0, 0, 0]

gen = ExtrudedSurfaceGenerator(surface, 
                               lconstant = 0.05,
                               lextrude = 0.5, 
                               nextrude = 20, 
                               dltarget = 0.005, 
                               dstarget = 0.1,
                               rho = 1.0,
                               flags = flags,
                               nNodePerh = 1.01,
                               SPH = False)

eos = GammaLawGasMKS(5.0/3.0, 1.0)
nodes = makeFluidNodeList("test", eos, 
                          hminratio = 1.0e-10)

distributeNodes3d((nodes, gen))

H = nodes.Hfield()
Hinv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    Hinv[i] = H[i].Inverse()

siloPointmeshDump("test_void", 
                  fields = [nodes.mass(), nodes.massDensity(), nodes.Hfield(), Hinv])
writePolyhedronOBJ(surface, "test_void_surface.obj")
