from Spheral3d import *
from RecursivePrimitiveRefinementGenerator import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from SpheralTestUtilities import *
from SpheralVisitDump import *


commandLine(nPerh   = 2.01,
            hmin    = 1e-5,
            hmax    = 1e6,
            rmin    = 0.0,
            rmax    = 2.0,
            scaler  = 0.6,
            n      = 2000)

class densityProf:
    def __init__(self,scaleR):
        self.R = scaleR
        return
    def __call__(self,r):
        return 1.0e6*exp(-r/self.R)
#return (2.2-r)/self.R

class rejecter(object):
    def __init__(self,point,radius):
        self.point=point
        self.radius=radius
    def accept(self,x,y,z):
        dx = self.point[0] - x
        dy = self.point[1] - y
        dz = self.point[2] - z
        if (sqrt(dx*dx+dy*dy+dz*dz)<self.radius):
            return True
        else:
            return False

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
reject = rejecter([2,0,0],0.5)
generator = RPRPSGenerator3d(n,rhoProfile,rmin,rmax,
                                                 nNodePerh = nPerh,rejecter=None)

nodes.numInternalNodes = generator.localNumNodes()

distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in range(nodes.numNodes):
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

procs = mpi.procs
rank = mpi.rank
serialData = []
i,j = 0,0
for i in range(procs):
    if rank == i:
        for j in range(nodes.numInternalNodes):
            serialData.append([nodes.positions()[j],3.0/(nodes.Hfield()[j].Trace()),nodes.mass()[j],nodes.massDensity()[j],nodes.specificThermalEnergy()[j]])
serialData = mpi.reduce(serialData,mpi.SUM)
if rank == 0:
    f = open("serialIcoDump.ascii",'w')
    for i in range(len(serialData)):
        f.write("{0} {1} {2} {3} {4}\n".format(i,sqrt(pow(serialData[i][0][0],2.0)+pow(serialData[i][0][1],2.0)+pow(serialData[i][0][2],2.0)),serialData[i][1],serialData[i][2],serialData[i][3]))
    f.close()

