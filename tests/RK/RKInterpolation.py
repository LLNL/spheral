#-------------------------------------------------------------------------------
# Manufactured diffusion test
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
title("Reproducing kernel test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial parameters
    dimension = 1,
    nx = 32,
    ny = 32,
    nz = 32,
    correctionOrder = LinearOrder,
    
    # Manufactured parameters
    rho0 = 2.5e-7, # g / cm^3
    
    # Interpolation kernel choice
    nPerh = 4.01,
    
    # Gamma-law gas EOS
    gamma = 5.0/3.0,
    mu = 1.0,

    # Data dir
    dataDirBase = "dumps-RKInterpolation"
)

if dimension == 1:
    from Spheral1d import *
elif dimension == 2:
    from Spheral2d import *
else:
    from Spheral3d import *

#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------

# Set limits to be two wavelengths so boundary conditions work
x0 = y0 = z0 = -1.0
x1 = y1 = z1 = 1.0
lims = [[x0, x1], [y0, y1], [z0, z1]]

# Get point spacing
units = MKS()
delta = [(x1-x0)/nx]
if dimension > 1:
    delta.append((y1 -y0)/ny)
if dimension > 2:
    delta.append((z1 - z0)/nz)
deltaMax = max(delta)
deltaMin = min(delta)
hmin = deltaMin * nPerh
hmax = deltaMax * nPerh
h0 = 0.5 * (hmin + hmax)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "dim={}".format(dimension),
                       "nx={}".format(nx))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "RKInterpolation")

import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh,
                          kernelExtent = kernelExtent)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Seed the nodes
#-------------------------------------------------------------------------------
if dimension == 1:
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (x0, x1))],
                             nPerh = nPerh)
elif dimension == 2:
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(distributionType="lattice",
                                           nRadial = nx,
                                           nTheta = ny,
                                           xmin = (x0, y0),
                                           xmax = (x1, y1),
                                           rho = rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
    distributeNodes2d((nodes, generator))
else:
    from GenerateNodeDistribution3d import *
    generator = GenerateNodeDistribution3d(distributionType="lattice",
                                           n1 = nx,
                                           n2 = ny,
                                           n3 = nz,
                                           xmin = (x0, y0, z0),
                                           xmax = (x1, y1, z1),
                                           rho=rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes, generator))
    
output("nodes.numNodes")
numLocal = nodes.numInternalNodes
output("numLocal")

#-------------------------------------------------------------------------------
# Build the RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(dataBase = db,
                   W = WT,
                   correctionOrder = correctionOrder,
                   volumeType = CRKMassOverDensity)
packages = [rk]

#-------------------------------------------------------------------------------
# Run the startup stuff 
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
control = SpheralController(integrator, WT)

#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
position = db.fluidPosition
