from Spheral3d import *
from GenerateNodeDistribution3d import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralVisitDump import *


commandLine(nPerh   = 2.01,
            hmin    = 1e-5,
            hmax    = 1e6,
            rmin    = 0.2,
            rmax    = 2.0,
            zmax    = 0.5,
            scaler  = 0.5,
            scalez  = 0.1,
            nr      = 100)

class diskDensity:
    def __init__(self,scaleR,scaleH):
        self.H = scaleH
        self.R = scaleR
        return
    def __call__(self,r,z):
        return exp(-r/self.R)*exp(-abs(z)/self.H)


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
densityProfile = diskDensity(scaler,scalez)

generator = GenerateIdealDiskMatchingProfile3d(nr,densityProfile,rmin,rmax,zmax,
                                               nNodePerh = nPerh)
nodes.numInternalNodes = generator.localNumNodes()

distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = SpheralVisitDump(baseFileName = "disk_test",
                           listOfFields = [nodes.massDensity(),
                                           nodes.mass(),
                                           nodes.velocity(),
                                           nodes.specificThermalEnergy(),
                                           Hfield,
                                           HfieldInv],
                           )
vizfile.dump(0.0, 0)
