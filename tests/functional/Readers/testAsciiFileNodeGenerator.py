from Spheral3d import *
from AsciiFileNodeGenerator import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralVisitDump import *

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
gamma = 1.4
mu = 2.0
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
hmin = 1e-5
hmax = 1e6
nPerh = 1.51
nodes = makeFluidNodeList("nodes", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh,
                          xmin = Vector.one * -1e20,
                          xmax = Vector.one * 1e20)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
generator = AsciiFileNodeGenerator3D(filename = "ic.sdf.ascii",
                                materialName = "Default",
                                nNodePerh = nPerh)
nodes.numInternalNodes = generator.localNumNodes()
vel = nodes.velocity()
eps = nodes.specificThermalEnergy()
abund = []

for i in range(nodes.numInternalNodes):
    vel[i].x = generator.vx[i]
    vel[i].y = generator.vy[i]
    vel[i].z = generator.vz[i]
    eps[i] = generator.eps[i]

distributeNodes((nodes, generator),)

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in range(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = SpheralVisitDump(baseFileName = "Ascii_file_test",
                           listOfFields = [nodes.massDensity(),
                                           nodes.mass(),
                                           nodes.velocity(),
                                           nodes.specificThermalEnergy(),
                                           Hfield,
                                           HfieldInv],
                           )
vizfile.dump(0.0, 0)
