from Spheral3d import *
from GenerateRatioSphere import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from siloPointmeshDump import *

commandLine(hmin    = 1e-5,
            hmax    = 1e6,
            rhoscale = 0.5,
            drCenter = 0.01,
            drRatio = 1.2,
            rmin    = 0.0,
            rmax    = 2.0,
            thetamin = 0.0,
            thetamax = 1.0*pi,
            phi = 0.75*pi,
            ntheta = 100,
            xcenter = 0.5,
            ycenter = 0.5,
            zcenter = 0.5,
            distributionType = "constantDTheta",
            nPerh   = 2.01,
            nr      = 100,
            SPH     = False)

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
                          topGridCellSize = 100,
                          xmin = Vector.one * -100.0,
                          xmax = Vector.one *  100.0)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
def rhoprofile(r):
    return exp(-r*r/(rhoscale*rhoscale))

generator = GenerateRatioSphere3d(drCenter, drRatio,
                                  rho = rhoprofile,
                                  rmin = rmin,
                                  rmax = rmax,
                                  thetamin = thetamin,
                                  thetamax = thetamax,
                                  phi = phi,
                                  ntheta = ntheta,
                                  center = (xcenter, ycenter, zcenter),
                                  distributionType = distributionType,
                                  nNodePerh = nPerh,
                                  SPH = SPH)
distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = SymTensor(Hfield[i].Inverse())
vizfile = siloPointmeshDump(baseName = "ratio_sphere_test_" + distributionType,
                            baseDirectory = "ratio_sphere_test_" + distributionType,
                            fields = [nodes.massDensity(),
                                      nodes.mass(),
                                      nodes.velocity(),
                                      nodes.specificThermalEnergy(),
                                      Hfield,
                                      HfieldInv],
                            )
