from Spheral3d import *
from GenerateEqualMassSheets3d import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralVisitDump import *


commandLine(nPerh   = 2.01,
            hmin    = 1e-5,
            hmax    = 1e6,
            rmin    = 0.0,
            rmax    = 2.0,
            scaler  = 0.6,
            nr      = 20,
            seed    = "ico")

class densityProf:
    def __init__(self,scaleR):
        self.R = scaleR
        return
    def __call__(self,r):
        return 1.0e6*exp(-r/self.R)
#return (2.2-r)/self.R


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
rhoProfile = densityProf(scaler)

generator = GenerateEqualMassSheets3d(nr,nr,rhoProfile,Vector.one * rmin,Vector.one * rmax,
                                      nNodePerh = nPerh,
                                      rhoMin = rhoProfile(rmax))

nodes.numInternalNodes = generator.localNumNodes()

distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = SpheralVisitDump(baseFileName = "icosahedron_test",
                           listOfFields = [nodes.massDensity(),
                                           nodes.mass(),
                                           nodes.velocity(),
                                           nodes.specificThermalEnergy(),
                                           Hfield,
                                           HfieldInv],
                           )
vizfile.dump(0.0, 0)

