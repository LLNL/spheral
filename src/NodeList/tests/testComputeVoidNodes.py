#-------------------------------------------------------------------------------
# Create some simple 3D node distrbibutions and test out creating void nodes
# on 'em.
#-------------------------------------------------------------------------------
from math import *
from Spheral3d import *
from GenerateNodeDistribution3d import *
from VoronoiDistributeNodes import distributeNodes3d
from generateMesh import *
from siloMeshDump import *

gamma = 5.0/3.0
mu = 1.0
r = 1.0
z0, z1, z2 = -1.0, 0.0, 1.0
nr = 20
nz = 10

eos = GammaLawGasMKS(gamma, mu)
nodes1 = makeFluidNodeList("nodes 1", eos)
nodes2 = makeFluidNodeList("nodes 2", eos)
nodes3 = makeFluidNodeList("nodes 3", eos)
generator1 = GenerateNodeDistribution3d(nr, nz, 0, 1.0, "cylindrical",
                                        rmin = 0.0,
                                        rmax = r,
                                        thetamin = 0.0,
                                        thetamax = 2.0*pi,
                                        zmin = z0,
                                        zmax = z1,
                                        nNodePerh = 2.01,
                                        SPH = True)
generator2 = GenerateNodeDistribution3d(nr, nz, 0, 1.0, "cylindrical",
                                        rmin = 0.0,
                                        rmax = r,
                                        thetamin = 0.0,
                                        thetamax = 2.0*pi,
                                        zmin = z1,
                                        zmax = z2,
                                        nNodePerh = 2.01,
                                        SPH = True)
distributeNodes3d((nodes1, generator1),
                  (nodes2, generator2))

mesh, void = generatePolyhedralMesh([nodes1, nodes2, nodes3],
                                    xmin = Vector(-2, -2, -2),
                                    xmax = Vector(2, 2, 2),
                                    generateVoid = False,
                                    generateParallelConnectivity = False,
                                    removeBoundaryZones = True)

nodeLists = vector_of_NodeList()
for nodes in [nodes1, nodes2, nodes3]:
    nodeLists.append(nodes)
    nodes.neighbor().updateNodes()
W = TableKernel(BSplineKernel(), 1000)
mom0 = ScalarFieldList(CopyFields)
mom1 = VectorFieldList(CopyFields)
zerothAndFirstNodalMoments(nodeLists, W, True, mom0, mom1)

siloMeshDump("voidNodeTest", mesh,
             nodeLists = [nodes1, nodes2, nodes3],
             scalarFields = list(mom0),
             vectorFields = list(mom1))
