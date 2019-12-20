#ATS:test(SELF, "--dimension 1 --correctionOrder ZerothOrder --funcType constant", label="RK interpolation - 1D zeroth")
#ATS:test(SELF, "--dimension 2 --correctionOrder ZerothOrder --funcType constant --numToCheck 10", label="RK interpolation - 2D zeroth")
#ATS:test(SELF, "--dimension 3 --correctionOrder ZerothOrder --funcType constant --numToCheck 10", label="RK interpolation - 3D zeroth")

#ATS:test(SELF, "--dimension 2 --correctionOrder LinearOrder --funcType linear", label="RK interpolation - 2D linear")
#ATS:test(SELF, "--dimension 3 --correctionOrder LinearOrder --funcType linear --numToCheck 10", label="RK interpolation - 3D linear")
#ATS:test(SELF, "--dimension 1 --correctionOrder LinearOrder --funcType linear --numToCheck 10", label="RK interpolation - 1D linear")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True", label="RK interpolation - 1D quadratic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True --numToCheck 10", label="RK interpolation - 2D quadratic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True --numToCheck 10", label="RK interpolation - 3D quadratic")

#ATS:test(SELF, "--dimension 1 --correctionOrder CubicOrder --funcType cubic --testHessian True", label="RK interpolation - 1D cubic")
#ATS:test(SELF, "--dimension 2 --correctionOrder CubicOrder --funcType cubic --testHessian True --numToCheck 10", label="RK interpolation - 2D cubic")
#ATS:test(SELF, "--dimension 3 --correctionOrder CubicOrder --funcType cubic --testHessian True --numToCheck 10", label="RK interpolation - 3D cubic")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True", label="RK interpolation - 1D quartic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True --numToCheck 10", label="RK interpolation - 2D quartic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True --numToCheck 10", label="RK interpolation - 3D quartic")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True", label="RK interpolation - 1D quintic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True --numToCheck 10", label="RK interpolation - 2D quintic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True --numToCheck 10", label="RK interpolation - 3D quintic")

#ATS:test(SELF, "--dimension 1 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True", label="RK interpolation - 1D sextic")
#ATS:test(SELF, "--dimension 2 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True --numToCheck 10", label="RK interpolation - 2D sextic")
#ATS:test(SELF, "--dimension 3 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True --numToCheck 10", label="RK interpolation - 3D sextic")

#ATS:test(SELF, "--dimension 1 --correctionOrder SepticOrder --funcType septic --nx 12 --nPerh 8.01 --tolerance 1.e-8 --testHessian True", label="RK interpolation - 1D septic")
#ATS:test(SELF, "--dimension 2 --correctionOrder SepticOrder --funcType septic --nPerh 8.01 --tolerance 1.e-8 --testHessian True --numToCheck 10", label="RK interpolation - 2D septic")
#ATS:test(SELF, "--dimension 3 --correctionOrder SepticOrder --funcType septic --nPerh 8.01 --tolerance 1.e-8 --testHessian True --numToCheck 10", label="RK interpolation - 3D septic")

#-------------------------------------------------------------------------------
# Test of interpolation for reproducing kernels
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
import time
from matplotlib import pyplot as plt
title("Reproducing kernel test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial parameters
    dimension = 1,
    nx = 10,
    ny = 10,
    nz = 10,
    x0 = -2.0,
    x1 = 2.0,
    y0 = -2.0,
    y1 = 2.0,
    z0 = -2.0,
    z1 = 2.0,

    # RK options
    correctionOrder = LinearOrder,
    volumeType = RKMassOverDensity,

    # Testing options
    randomizeNodes = True,
    ranfrac = 0.2,
    testHessian = False,
    useOldKernel = False, # Test using old kernel
    useBaseKernel = False, # Test using standard SPH kernel
    checkAgainstOld = False,
    printErrors = False,
    quitAfterTiming = False,
    computeCorrectionsDirectly = False,
    numToCheck = -1, # -1 for all nodes, positive int for a few random nodes
    
    # Manufactured parameters
    funcType = "linear",
    
    # Interpolation kernel choice
    nPerh = 4.01,
    hminmult = 1.e-3,
    hmaxmult = 1.e3,

    # Options for checking
    checkConditions = True,
    
    # Material parameters
    rho0 = 2.5e-7,
    gamma = 5.0/3.0,
    mu = 1.0,

    # Data dir
    dataDirBase = "dumps-RKInterpolation",

    # Error
    tolerance = 1.e-12,
    
    # Plotting
    plot = False,
)

if useOldKernel and useBaseKernel:
    raise ValueError, "cannot use both old and base kernel"
if useOldKernel and correctionOrder > QuadraticOrder:
    raise ValueError, "correction order must be quadratic to use old kernel"

if nPerh < int(correctionOrder):
    print "nPerh is not large enough for correction order: {} < {}".format(nPerh, int(correctionOrder))
    
if mpi.procs > 1:
    raise ValueError, "parallel node generation not working"
    
#-------------------------------------------------------------------------------
# Choose correct corrections
#-------------------------------------------------------------------------------
if dimension == 1:
    from Spheral1d import *
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections1d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections1d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections1d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections1d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections1d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections1d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections1d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections1d7
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)
elif dimension == 2:
    from Spheral2d import *
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections2d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections2d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections2d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections2d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections2d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections2d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections2d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections2d7
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)
else:
    from Spheral3d import *
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections3d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections3d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections3d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections3d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections3d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections3d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections3d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections3d7
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)

RKUtilities = makeRKUtilities(correctionOrder)

#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------
# Set limits
lims = [[x0, x1], [y0, y1], [z0, z1]]

# Get point spacing
units = MKS()
delta = [(x1-x0)/nx]
if dimension > 1:
    delta.append((y1 -y0)/ny)
if dimension > 2:
    delta.append((z1 - z0)/nz)
deltaMin = min(delta)
deltaMax = max(delta)
hmin = deltaMin * nPerh * hminmult
hmax = deltaMax * nPerh * hmaxmult

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "correctionOrder={}".format(correctionOrder),
                       "funcType={}".format(funcType),
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
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")
output("dataBase.numNodeLists")
output("dataBase.numFluidNodeLists")

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
# Randomize nodes
#-------------------------------------------------------------------------------
import random
seed = 2
rangen = random.Random()
rangen.seed(seed)

if randomizeNodes:
    dx = (x1 - x0)/nx
    dy = (y1 - y0)/ny
    dz = (z1 - z0)/nz
    pos = nodes.positions()
    for i in xrange(nodes.numInternalNodes):
        if dimension == 1:
            pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
        elif dimension == 2:
            pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * rangen.uniform(-1.0, 1.0)
        elif dimension == 3:
            pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * rangen.uniform(-1.0, 1.0)
            pos[i].z += ranfrac * dz * rangen.uniform(-1.0, 1.0)
            
#-------------------------------------------------------------------------------
# Iterate h
#-------------------------------------------------------------------------------
bounds = vector_of_Boundary()
method = SPHSmoothingScale()
iterateIdealH(dataBase,
              bounds,
              WT,
              method,
              100, # max h iterations
              1.e-4) # h tolerance
dataBase.updateConnectivityMap(True)

#-------------------------------------------------------------------------------
# Get interpolant
#-------------------------------------------------------------------------------
if funcType == "constant":
    def func(x):
        return 2.0
    def dfunc(x):
        return [0.0 for v in x]
    def ddfunc(x):
        return [[0.0 for v in x] for v in x]
elif funcType == "linear":
    if dimension == 1:
        def func(x):
            return -10 - 5*x[0]
        def dfunc(x):
            return [-5]
        def ddfunc(x):
            return [[0]]
    elif dimension == 2:
        def func(x):
            return 5 + 8*x[0] - 4*x[1]
        def dfunc(x):
            return [8, -4]
        def ddfunc(x):
            return [[0, 0], [0, 0]]
    else:
        def func(x):
            return -5 + 2*x[0] - 4*x[1] - 9*x[2]
        def dfunc(x):
            return [2, -4, -9]
        def ddfunc(x):
            return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
elif funcType == "quadratic":
    if dimension == 1:
        def func(x):
            return 9 - 9*x[0] - 9*np.power(x[0],2)
        def dfunc(x):
            return [-9 - 18*x[0]]
        def ddfunc(x):
            return [[-18]]
    elif dimension == 2:
        def func(x):
            return 4 - 7*x[0] - np.power(x[0],2) - 6*x[1] + 9*x[0]*x[1] + 10*np.power(x[1],2)
        def dfunc(x):
            return [-7 - 2*x[0] + 9*x[1], -6 + 9*x[0] + 20*x[1]]
        def ddfunc(x):
            return [[-2, 9], [9, 20]]
    else:
        def func(x):
            return 10 - 6*x[0] - np.power(x[0],2) - 6*x[1] + 4*x[0]*x[1] + np.power(x[1],2) + 8*x[2] + 6*x[0]*x[2] - 7*x[1]*x[2] + 3*np.power(x[2],2)
        def dfunc(x):
            return [-6 - 2*x[0] + 4*x[1] + 6*x[2], -6 + 4*x[0] + 2*x[1] - 7*x[2], 8 + 6*x[0] - 7*x[1] + 6*x[2]]
        def ddfunc(x):
            return [[-2, 4, 6], [4, 2, -7], [6, -7, 6]]
elif funcType == "cubic":
    if dimension == 1:
        def func(x):
            return 9 + 5*x[0] + 8*np.power(x[0],2) - 9*np.power(x[0],3)
        def dfunc(x):
            return [5 + 16*x[0] - 27*np.power(x[0],2)]
        def ddfunc(x):
            return [[16 - 54*x[0]]]
    elif dimension == 2:
        def func(x):
            return 1 + 10*x[0] - 4*np.power(x[0],2) + 8*np.power(x[0],3) + 10*x[1] - 2*x[0]*x[1] + 5*np.power(x[0],2)*x[1] + 6*np.power(x[1],2) + 8*x[0]*np.power(x[1],2) - 8*np.power(x[1],3)
        def dfunc(x):
            return [10 - 8*x[0] + 24*np.power(x[0],2) - 2*x[1] + 10*x[0]*x[1] + 8*np.power(x[1],2), 10 - 2*x[0] + 5*np.power(x[0],2) + 12*x[1] + 16*x[0]*x[1] - 24*np.power(x[1],2)]
        def ddfunc(x):
            return [[-8 + 48*x[0] + 10*x[1], -2 + 10*x[0] + 16*x[1]], [-2 + 10*x[0] + 16*x[1], 12 + 16*x[0] - 48*x[1]]]
    else:
        def func(x):
            return 6 + 8*x[0] + 5*np.power(x[0],2) + 2*np.power(x[0],3) - 5*x[1] - 8*x[0]*x[1] - np.power(x[0],2)*x[1] - 10*np.power(x[1],2) - 9*x[0]*np.power(x[1],2) - 4*np.power(x[1],3) - 8*x[2] + 9*x[0]*x[2] - 10*np.power(x[0],2)*x[2] - 7*x[1]*x[2] - 7*x[0]*x[1]*x[2] - 5*np.power(x[1],2)*x[2] + 2*np.power(x[2],2) + 5*x[0]*np.power(x[2],2) - 4*x[1]*np.power(x[2],2) + 10*np.power(x[2],3)
        def dfunc(x):
            return [8 + 10*x[0] + 6*np.power(x[0],2) - 8*x[1] - 2*x[0]*x[1] - 9*np.power(x[1],2) + 9*x[2] - 20*x[0]*x[2] - 7*x[1]*x[2] + 5*np.power(x[2],2), -5 - 8*x[0] - np.power(x[0],2) - 20*x[1] - 18*x[0]*x[1] - 12*np.power(x[1],2) - 7*x[2] - 7*x[0]*x[2] - 10*x[1]*x[2] - 4*np.power(x[2],2), -8 + 9*x[0] - 10*np.power(x[0],2) - 7*x[1] - 7*x[0]*x[1] - 5*np.power(x[1],2) + 4*x[2] + 10*x[0]*x[2] - 8*x[1]*x[2] + 30*np.power(x[2],2)]
        def ddfunc(x):
            return [[10 + 12*x[0] - 2*x[1] - 20*x[2], -8 - 2*x[0] - 18*x[1] - 7*x[2], 9 - 20*x[0] - 7*x[1] + 10*x[2]], [-8 - 2*x[0] - 18*x[1] - 7*x[2], -20 - 18*x[0] - 24*x[1] - 10*x[2], -7 - 7*x[0] - 10*x[1] - 8*x[2]], [9 - 20*x[0] - 7*x[1] + 10*x[2], -7 - 7*x[0] - 10*x[1] - 8*x[2], 4 + 10*x[0] - 8*x[1] + 60*x[2]]]
elif funcType == "quartic":
    if dimension == 1:
        def func(x):
            return 10 + 7*x[0] + 3*np.power(x[0],2) - 5*np.power(x[0],3) - 9*np.power(x[0],4)
        def dfunc(x):
            return [7 + 6*x[0] - 15*np.power(x[0],2) - 36*np.power(x[0],3)]
        def ddfunc(x):
            return [[6 - 30*x[0] - 108*np.power(x[0],2)]]
    elif dimension == 2:
        def func(x):
            return 4 + 9*x[0] - 2*np.power(x[0],2) + 3*np.power(x[0],3) + 2*np.power(x[0],4) - 9*x[1] + 7*x[0]*x[1] - 2*np.power(x[0],2)*x[1] - 3*np.power(x[0],3)*x[1] - 7*np.power(x[1],2) - 2*x[0]*np.power(x[1],2) - 2*np.power(x[0],2)*np.power(x[1],2) - 2*np.power(x[1],3) + 4*x[0]*np.power(x[1],3) + 4*np.power(x[1],4)
        def dfunc(x):
            return [9 - 4*x[0] + 9*np.power(x[0],2) + 8*np.power(x[0],3) + 7*x[1] - 4*x[0]*x[1] - 9*np.power(x[0],2)*x[1] - 2*np.power(x[1],2) - 4*x[0]*np.power(x[1],2) + 4*np.power(x[1],3), -9 + 7*x[0] - 2*np.power(x[0],2) - 3*np.power(x[0],3) - 14*x[1] - 4*x[0]*x[1] - 4*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) + 12*x[0]*np.power(x[1],2) + 16*np.power(x[1],3)]
        def ddfunc(x):
            return [[-4 + 18*x[0] + 24*np.power(x[0],2) - 4*x[1] - 18*x[0]*x[1] - 4*np.power(x[1],2), 7 - 4*x[0] - 9*np.power(x[0],2) - 4*x[1] - 8*x[0]*x[1] + 12*np.power(x[1],2)], [7 - 4*x[0] - 9*np.power(x[0],2) - 4*x[1] - 8*x[0]*x[1] + 12*np.power(x[1],2), -14 - 4*x[0] - 4*np.power(x[0],2) - 12*x[1] + 24*x[0]*x[1] + 48*np.power(x[1],2)]]
    else:
        def func(x):
            return 9 - x[0] - np.power(x[0],2) + 8*np.power(x[0],3) - 7*np.power(x[0],4) - 7*x[1] - 3*x[0]*x[1] - 3*np.power(x[0],2)*x[1] + 6*np.power(x[0],3)*x[1] + 5*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 7*np.power(x[0],2)*np.power(x[1],2) + np.power(x[1],3) - 2*x[0]*np.power(x[1],3) + 4*np.power(x[1],4) - x[2] + 8*x[0]*x[2] - 8*np.power(x[0],2)*x[2] + 4*np.power(x[0],3)*x[2] - 6*x[1]*x[2] + 2*x[0]*x[1]*x[2] - 2*np.power(x[0],2)*x[1]*x[2] + 7*np.power(x[1],2)*x[2] - 5*x[0]*np.power(x[1],2)*x[2] - 9*np.power(x[1],3)*x[2] - np.power(x[2],2) - 4*x[0]*np.power(x[2],2) - 2*np.power(x[0],2)*np.power(x[2],2) + 2*x[1]*np.power(x[2],2) + 2*x[0]*x[1]*np.power(x[2],2) + 8*np.power(x[1],2)*np.power(x[2],2) - 7*np.power(x[2],3) - 10*x[0]*np.power(x[2],3) + 2*x[1]*np.power(x[2],3) + 7*np.power(x[2],4)
        def dfunc(x):
            return [-1 - 2*x[0] + 24*np.power(x[0],2) - 28*np.power(x[0],3) - 3*x[1] - 6*x[0]*x[1] + 18*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 14*x[0]*np.power(x[1],2) - 2*np.power(x[1],3) + 8*x[2] - 16*x[0]*x[2] + 12*np.power(x[0],2)*x[2] + 2*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 5*np.power(x[1],2)*x[2] - 4*np.power(x[2],2) - 4*x[0]*np.power(x[2],2) + 2*x[1]*np.power(x[2],2) - 10*np.power(x[2],3), -7 - 3*x[0] - 3*np.power(x[0],2) + 6*np.power(x[0],3) + 10*x[1] - 12*x[0]*x[1] - 14*np.power(x[0],2)*x[1] + 3*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) + 16*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] - 2*np.power(x[0],2)*x[2] + 14*x[1]*x[2] - 10*x[0]*x[1]*x[2] - 27*np.power(x[1],2)*x[2] + 2*np.power(x[2],2) + 2*x[0]*np.power(x[2],2) + 16*x[1]*np.power(x[2],2) + 2*np.power(x[2],3), -1 + 8*x[0] - 8*np.power(x[0],2) + 4*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] - 2*np.power(x[0],2)*x[1] + 7*np.power(x[1],2) - 5*x[0]*np.power(x[1],2) - 9*np.power(x[1],3) - 2*x[2] - 8*x[0]*x[2] - 4*np.power(x[0],2)*x[2] + 4*x[1]*x[2] + 4*x[0]*x[1]*x[2] + 16*np.power(x[1],2)*x[2] - 21*np.power(x[2],2) - 30*x[0]*np.power(x[2],2) + 6*x[1]*np.power(x[2],2) + 28*np.power(x[2],3)]
        def ddfunc(x):
            return [[-2 + 48*x[0] - 84*np.power(x[0],2) - 6*x[1] + 36*x[0]*x[1] - 14*np.power(x[1],2) - 16*x[2] + 24*x[0]*x[2] - 4*x[1]*x[2] - 4*np.power(x[2],2), -3 - 6*x[0] + 18*np.power(x[0],2) - 12*x[1] - 28*x[0]*x[1] - 6*np.power(x[1],2) + 2*x[2] - 4*x[0]*x[2] - 10*x[1]*x[2] + 2*np.power(x[2],2), 8 - 16*x[0] + 12*np.power(x[0],2) + 2*x[1] - 4*x[0]*x[1] - 5*np.power(x[1],2) - 8*x[2] - 8*x[0]*x[2] + 4*x[1]*x[2] - 30*np.power(x[2],2)], [-3 - 6*x[0] + 18*np.power(x[0],2) - 12*x[1] - 28*x[0]*x[1] - 6*np.power(x[1],2) + 2*x[2] - 4*x[0]*x[2] - 10*x[1]*x[2] + 2*np.power(x[2],2), 10 - 12*x[0] - 14*np.power(x[0],2) + 6*x[1] - 12*x[0]*x[1] + 48*np.power(x[1],2) + 14*x[2] - 10*x[0]*x[2] - 54*x[1]*x[2] + 16*np.power(x[2],2), -6 + 2*x[0] - 2*np.power(x[0],2) + 14*x[1] - 10*x[0]*x[1] - 27*np.power(x[1],2) + 4*x[2] + 4*x[0]*x[2] + 32*x[1]*x[2] + 6*np.power(x[2],2)], [8 - 16*x[0] + 12*np.power(x[0],2) + 2*x[1] - 4*x[0]*x[1] - 5*np.power(x[1],2) - 8*x[2] - 8*x[0]*x[2] + 4*x[1]*x[2] - 30*np.power(x[2],2), -6 + 2*x[0] - 2*np.power(x[0],2) + 14*x[1] - 10*x[0]*x[1] - 27*np.power(x[1],2) + 4*x[2] + 4*x[0]*x[2] + 32*x[1]*x[2] + 6*np.power(x[2],2), -2 - 8*x[0] - 4*np.power(x[0],2) + 4*x[1] + 4*x[0]*x[1] + 16*np.power(x[1],2) - 42*x[2] - 60*x[0]*x[2] + 12*x[1]*x[2] + 84*np.power(x[2],2)]]
elif funcType == "quintic":
    if dimension == 1:
        def func(x):
            return -3 - 10*x[0] + np.power(x[0],2) + 8*np.power(x[0],3) - 6*np.power(x[0],4) + 5*np.power(x[0],5)
        def dfunc(x):
            return [-10 + 2*x[0] + 24*np.power(x[0],2) - 24*np.power(x[0],3) + 25*np.power(x[0],4)]
        def ddfunc(x):
            return [[2 + 48*x[0] - 72*np.power(x[0],2) + 100*np.power(x[0],3)]]
    elif dimension == 2:
        def func(x):
            return 4 + 6*x[0] - 7*np.power(x[0],2) - 2*np.power(x[0],3) + np.power(x[0],4) + 4*np.power(x[0],5) + 2*x[1] - 7*x[0]*x[1] + 5*np.power(x[0],2)*x[1] - np.power(x[0],3)*x[1] - 6*np.power(x[0],4)*x[1] + 10*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 4*np.power(x[0],2)*np.power(x[1],2) + np.power(x[0],3)*np.power(x[1],2) - 8*np.power(x[1],3) - 2*x[0]*np.power(x[1],3) - 10*np.power(x[0],2)*np.power(x[1],3) + 2*np.power(x[1],4) - 7*x[0]*np.power(x[1],4) - 8*np.power(x[1],5)
        def dfunc(x):
            return [6 - 14*x[0] - 6*np.power(x[0],2) + 4*np.power(x[0],3) + 20*np.power(x[0],4) - 7*x[1] + 10*x[0]*x[1] - 3*np.power(x[0],2)*x[1] - 24*np.power(x[0],3)*x[1] - 6*np.power(x[1],2) - 8*x[0]*np.power(x[1],2) + 3*np.power(x[0],2)*np.power(x[1],2) - 2*np.power(x[1],3) - 20*x[0]*np.power(x[1],3) - 7*np.power(x[1],4), 2 - 7*x[0] + 5*np.power(x[0],2) - np.power(x[0],3) - 6*np.power(x[0],4) + 20*x[1] - 12*x[0]*x[1] - 8*np.power(x[0],2)*x[1] + 2*np.power(x[0],3)*x[1] - 24*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 30*np.power(x[0],2)*np.power(x[1],2) + 8*np.power(x[1],3) - 28*x[0]*np.power(x[1],3) - 40*np.power(x[1],4)]
        def ddfunc(x):
            return [[-14 - 12*x[0] + 12*np.power(x[0],2) + 80*np.power(x[0],3) + 10*x[1] - 6*x[0]*x[1] - 72*np.power(x[0],2)*x[1] - 8*np.power(x[1],2) + 6*x[0]*np.power(x[1],2) - 20*np.power(x[1],3), -7 + 10*x[0] - 3*np.power(x[0],2) - 24*np.power(x[0],3) - 12*x[1] - 16*x[0]*x[1] + 6*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 60*x[0]*np.power(x[1],2) - 28*np.power(x[1],3)], [-7 + 10*x[0] - 3*np.power(x[0],2) - 24*np.power(x[0],3) - 12*x[1] - 16*x[0]*x[1] + 6*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 60*x[0]*np.power(x[1],2) - 28*np.power(x[1],3), 20 - 12*x[0] - 8*np.power(x[0],2) + 2*np.power(x[0],3) - 48*x[1] - 12*x[0]*x[1] - 60*np.power(x[0],2)*x[1] + 24*np.power(x[1],2) - 84*x[0]*np.power(x[1],2) - 160*np.power(x[1],3)]]
    else:
        def func(x):
            return -8 + 10*x[0] + 4*np.power(x[0],2) + 4*np.power(x[0],3) + 6*np.power(x[0],4) + 9*np.power(x[0],5) + 9*x[1] - 10*x[0]*x[1] + 2*np.power(x[0],2)*x[1] + np.power(x[0],3)*x[1] + 6*np.power(x[0],4)*x[1] - 8*np.power(x[1],2) - 4*x[0]*np.power(x[1],2) + np.power(x[0],2)*np.power(x[1],2) + 2*np.power(x[0],3)*np.power(x[1],2) + 10*np.power(x[1],3) - 10*x[0]*np.power(x[1],3) - 7*np.power(x[0],2)*np.power(x[1],3) - 10*np.power(x[1],4) + 10*x[0]*np.power(x[1],4) - 2*np.power(x[1],5) + 10*x[2] - 4*x[0]*x[2] - 5*np.power(x[0],2)*x[2] - 5*np.power(x[0],3)*x[2] + 2*np.power(x[0],4)*x[2] + 8*x[1]*x[2] - 6*x[0]*x[1]*x[2] + np.power(x[0],2)*x[1]*x[2] + 5*np.power(x[0],3)*x[1]*x[2] - 7*np.power(x[1],2)*x[2] - x[0]*np.power(x[1],2)*x[2] - 8*np.power(x[0],2)*np.power(x[1],2)*x[2] - 5*np.power(x[1],3)*x[2] + 8*x[0]*np.power(x[1],3)*x[2] + 10*np.power(x[1],4)*x[2] - 10*np.power(x[2],2) - 3*x[0]*np.power(x[2],2) + 6*np.power(x[0],2)*np.power(x[2],2) - 7*np.power(x[0],3)*np.power(x[2],2) - 2*x[1]*np.power(x[2],2) - 6*x[0]*x[1]*np.power(x[2],2) + 8*np.power(x[0],2)*x[1]*np.power(x[2],2) + 6*np.power(x[1],2)*np.power(x[2],2) - 7*x[0]*np.power(x[1],2)*np.power(x[2],2) + 2*np.power(x[1],3)*np.power(x[2],2) - 5*np.power(x[2],3) - 4*x[0]*np.power(x[2],3) - 7*np.power(x[0],2)*np.power(x[2],3) - 4*x[1]*np.power(x[2],3) + x[0]*x[1]*np.power(x[2],3) - 6*np.power(x[1],2)*np.power(x[2],3) - 3*np.power(x[2],4) - 10*x[0]*np.power(x[2],4) + 7*x[1]*np.power(x[2],4) - 8*np.power(x[2],5)
        def dfunc(x):
            return [10 + 8*x[0] + 12*np.power(x[0],2) + 24*np.power(x[0],3) + 45*np.power(x[0],4) - 10*x[1] + 4*x[0]*x[1] + 3*np.power(x[0],2)*x[1] + 24*np.power(x[0],3)*x[1] - 4*np.power(x[1],2) + 2*x[0]*np.power(x[1],2) + 6*np.power(x[0],2)*np.power(x[1],2) - 10*np.power(x[1],3) - 14*x[0]*np.power(x[1],3) + 10*np.power(x[1],4) - 4*x[2] - 10*x[0]*x[2] - 15*np.power(x[0],2)*x[2] + 8*np.power(x[0],3)*x[2] - 6*x[1]*x[2] + 2*x[0]*x[1]*x[2] + 15*np.power(x[0],2)*x[1]*x[2] - np.power(x[1],2)*x[2] - 16*x[0]*np.power(x[1],2)*x[2] + 8*np.power(x[1],3)*x[2] - 3*np.power(x[2],2) + 12*x[0]*np.power(x[2],2) - 21*np.power(x[0],2)*np.power(x[2],2) - 6*x[1]*np.power(x[2],2) + 16*x[0]*x[1]*np.power(x[2],2) - 7*np.power(x[1],2)*np.power(x[2],2) - 4*np.power(x[2],3) - 14*x[0]*np.power(x[2],3) + x[1]*np.power(x[2],3) - 10*np.power(x[2],4), 9 - 10*x[0] + 2*np.power(x[0],2) + np.power(x[0],3) + 6*np.power(x[0],4) - 16*x[1] - 8*x[0]*x[1] + 2*np.power(x[0],2)*x[1] + 4*np.power(x[0],3)*x[1] + 30*np.power(x[1],2) - 30*x[0]*np.power(x[1],2) - 21*np.power(x[0],2)*np.power(x[1],2) - 40*np.power(x[1],3) + 40*x[0]*np.power(x[1],3) - 10*np.power(x[1],4) + 8*x[2] - 6*x[0]*x[2] + np.power(x[0],2)*x[2] + 5*np.power(x[0],3)*x[2] - 14*x[1]*x[2] - 2*x[0]*x[1]*x[2] - 16*np.power(x[0],2)*x[1]*x[2] - 15*np.power(x[1],2)*x[2] + 24*x[0]*np.power(x[1],2)*x[2] + 40*np.power(x[1],3)*x[2] - 2*np.power(x[2],2) - 6*x[0]*np.power(x[2],2) + 8*np.power(x[0],2)*np.power(x[2],2) + 12*x[1]*np.power(x[2],2) - 14*x[0]*x[1]*np.power(x[2],2) + 6*np.power(x[1],2)*np.power(x[2],2) - 4*np.power(x[2],3) + x[0]*np.power(x[2],3) - 12*x[1]*np.power(x[2],3) + 7*np.power(x[2],4), 10 - 4*x[0] - 5*np.power(x[0],2) - 5*np.power(x[0],3) + 2*np.power(x[0],4) + 8*x[1] - 6*x[0]*x[1] + np.power(x[0],2)*x[1] + 5*np.power(x[0],3)*x[1] - 7*np.power(x[1],2) - x[0]*np.power(x[1],2) - 8*np.power(x[0],2)*np.power(x[1],2) - 5*np.power(x[1],3) + 8*x[0]*np.power(x[1],3) + 10*np.power(x[1],4) - 20*x[2] - 6*x[0]*x[2] + 12*np.power(x[0],2)*x[2] - 14*np.power(x[0],3)*x[2] - 4*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 16*np.power(x[0],2)*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 14*x[0]*np.power(x[1],2)*x[2] + 4*np.power(x[1],3)*x[2] - 15*np.power(x[2],2) - 12*x[0]*np.power(x[2],2) - 21*np.power(x[0],2)*np.power(x[2],2) - 12*x[1]*np.power(x[2],2) + 3*x[0]*x[1]*np.power(x[2],2) - 18*np.power(x[1],2)*np.power(x[2],2) - 12*np.power(x[2],3) - 40*x[0]*np.power(x[2],3) + 28*x[1]*np.power(x[2],3) - 40*np.power(x[2],4)]
        def ddfunc(x):
            return [[8 + 24*x[0] + 72*np.power(x[0],2) + 180*np.power(x[0],3) + 4*x[1] + 6*x[0]*x[1] + 72*np.power(x[0],2)*x[1] + 2*np.power(x[1],2) + 12*x[0]*np.power(x[1],2) - 14*np.power(x[1],3) - 10*x[2] - 30*x[0]*x[2] + 24*np.power(x[0],2)*x[2] + 2*x[1]*x[2] + 30*x[0]*x[1]*x[2] - 16*np.power(x[1],2)*x[2] + 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 16*x[1]*np.power(x[2],2) - 14*np.power(x[2],3), -10 + 4*x[0] + 3*np.power(x[0],2) + 24*np.power(x[0],3) - 8*x[1] + 4*x[0]*x[1] + 12*np.power(x[0],2)*x[1] - 30*np.power(x[1],2) - 42*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] + 15*np.power(x[0],2)*x[2] - 2*x[1]*x[2] - 32*x[0]*x[1]*x[2] + 24*np.power(x[1],2)*x[2] - 6*np.power(x[2],2) + 16*x[0]*np.power(x[2],2) - 14*x[1]*np.power(x[2],2) + np.power(x[2],3), -4 - 10*x[0] - 15*np.power(x[0],2) + 8*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] + 15*np.power(x[0],2)*x[1] - np.power(x[1],2) - 16*x[0]*np.power(x[1],2) + 8*np.power(x[1],3) - 6*x[2] + 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 12*x[1]*x[2] + 32*x[0]*x[1]*x[2] - 14*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 3*x[1]*np.power(x[2],2) - 40*np.power(x[2],3)], [-10 + 4*x[0] + 3*np.power(x[0],2) + 24*np.power(x[0],3) - 8*x[1] + 4*x[0]*x[1] + 12*np.power(x[0],2)*x[1] - 30*np.power(x[1],2) - 42*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] + 15*np.power(x[0],2)*x[2] - 2*x[1]*x[2] - 32*x[0]*x[1]*x[2] + 24*np.power(x[1],2)*x[2] - 6*np.power(x[2],2) + 16*x[0]*np.power(x[2],2) - 14*x[1]*np.power(x[2],2) + np.power(x[2],3), -16 - 8*x[0] + 2*np.power(x[0],2) + 4*np.power(x[0],3) + 60*x[1] - 60*x[0]*x[1] - 42*np.power(x[0],2)*x[1] - 120*np.power(x[1],2) + 120*x[0]*np.power(x[1],2) - 40*np.power(x[1],3) - 14*x[2] - 2*x[0]*x[2] - 16*np.power(x[0],2)*x[2] - 30*x[1]*x[2] + 48*x[0]*x[1]*x[2] + 120*np.power(x[1],2)*x[2] + 12*np.power(x[2],2) - 14*x[0]*np.power(x[2],2) + 12*x[1]*np.power(x[2],2) - 12*np.power(x[2],3), 8 - 6*x[0] + np.power(x[0],2) + 5*np.power(x[0],3) - 14*x[1] - 2*x[0]*x[1] - 16*np.power(x[0],2)*x[1] - 15*np.power(x[1],2) + 24*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 4*x[2] - 12*x[0]*x[2] + 16*np.power(x[0],2)*x[2] + 24*x[1]*x[2] - 28*x[0]*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) + 3*x[0]*np.power(x[2],2) - 36*x[1]*np.power(x[2],2) + 28*np.power(x[2],3)], [-4 - 10*x[0] - 15*np.power(x[0],2) + 8*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] + 15*np.power(x[0],2)*x[1] - np.power(x[1],2) - 16*x[0]*np.power(x[1],2) + 8*np.power(x[1],3) - 6*x[2] + 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 12*x[1]*x[2] + 32*x[0]*x[1]*x[2] - 14*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 3*x[1]*np.power(x[2],2) - 40*np.power(x[2],3), 8 - 6*x[0] + np.power(x[0],2) + 5*np.power(x[0],3) - 14*x[1] - 2*x[0]*x[1] - 16*np.power(x[0],2)*x[1] - 15*np.power(x[1],2) + 24*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 4*x[2] - 12*x[0]*x[2] + 16*np.power(x[0],2)*x[2] + 24*x[1]*x[2] - 28*x[0]*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) + 3*x[0]*np.power(x[2],2) - 36*x[1]*np.power(x[2],2) + 28*np.power(x[2],3), -20 - 6*x[0] + 12*np.power(x[0],2) - 14*np.power(x[0],3) - 4*x[1] - 12*x[0]*x[1] + 16*np.power(x[0],2)*x[1] + 12*np.power(x[1],2) - 14*x[0]*np.power(x[1],2) + 4*np.power(x[1],3) - 30*x[2] - 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 24*x[1]*x[2] + 6*x[0]*x[1]*x[2] - 36*np.power(x[1],2)*x[2] - 36*np.power(x[2],2) - 120*x[0]*np.power(x[2],2) + 84*x[1]*np.power(x[2],2) - 160*np.power(x[2],3)]]
elif funcType == "sextic":
    if dimension == 1:
        def func(x):
            return 5 - 3*x[0] - 8*(x[0]*x[0]) - 9*(x[0]*x[0]*x[0]) - 7*(x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])
        def dfunc(x):
            return [-3 - 16*x[0] - 27*(x[0]*x[0]) - 28*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0])]
        def ddfunc(x):
            return [[-16 - 54*x[0] - 84*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 270*(x[0]*x[0]*x[0]*x[0])]]
    elif dimension == 2:
        def func(x):
            return -1 + x[0] - x[0]*x[0] - 6*(x[0]*x[0]*x[0]) + 3*(x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] + 9*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - 6*(x[0]*x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 9*(x[1]*x[1]) + 6*x[0]*(x[1]*x[1]) + 6*(x[0]*x[0])*(x[1]*x[1]) + 5*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 10*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 9*(x[1]*x[1]*x[1]) - 7*x[0]*(x[1]*x[1]*x[1]) + 4*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]) + 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 7*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) - 2*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1]*x[1]
        def dfunc(x):
            return [1 - 2*x[0] - 18*(x[0]*x[0]) + 12*(x[0]*x[0]*x[0]) - 30*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*x[1] - 8*x[0]*x[1] - 18*(x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0])*x[1] - 45*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) + 12*x[0]*(x[1]*x[1]) + 15*(x[0]*x[0])*(x[1]*x[1]) + 40*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 7*(x[1]*x[1]*x[1]) + 8*x[0]*(x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 4*(x[1]*x[1]*x[1]*x[1]) + 14*x[0]*(x[1]*x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]*x[1]), -3 + 9*x[0] - 4*(x[0]*x[0]) - 6*(x[0]*x[0]*x[0]) - 4*(x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 18*x[1] + 12*x[0]*x[1] + 12*(x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0]*x[0])*x[1] + 27*(x[1]*x[1]) - 21*x[0]*(x[1]*x[1]) + 12*(x[0]*x[0])*(x[1]*x[1]) - 9*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 36*(x[1]*x[1]*x[1]) + 16*x[0]*(x[1]*x[1]*x[1]) + 28*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1])]
        def ddfunc(x):
            return [[-2 - 36*x[0] + 36*(x[0]*x[0]) - 120*(x[0]*x[0]*x[0]) + 270*(x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 36*x[0]*x[1] - 48*(x[0]*x[0])*x[1] - 180*(x[0]*x[0]*x[0])*x[1] + 12*(x[1]*x[1]) + 30*x[0]*(x[1]*x[1]) + 120*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]) + 14*(x[1]*x[1]*x[1]*x[1]), 9 - 8*x[0] - 18*(x[0]*x[0]) - 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) + 12*x[1] + 24*x[0]*x[1] + 30*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 21*(x[1]*x[1]) + 24*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 16*(x[1]*x[1]*x[1]) + 56*x[0]*(x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1])], [9 - 8*x[0] - 18*(x[0]*x[0]) - 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) + 12*x[1] + 24*x[0]*x[1] + 30*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 21*(x[1]*x[1]) + 24*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 16*(x[1]*x[1]*x[1]) + 56*x[0]*(x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1]), 18 + 12*x[0] + 12*(x[0]*x[0]) + 10*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) + 54*x[1] - 42*x[0]*x[1] + 24*(x[0]*x[0])*x[1] - 18*(x[0]*x[0]*x[0])*x[1] + 108*(x[1]*x[1]) + 48*x[0]*(x[1]*x[1]) + 84*(x[0]*x[0])*(x[1]*x[1]) - 40*(x[1]*x[1]*x[1]) - 80*x[0]*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1])]]
    else:
        def func(x):
            return 7 - 7*x[0] - 2*(x[0]*x[0]) - 3*(x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0] - 10*(x[0]*x[0]*x[0]*x[0]*x[0]) - 7*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] - 5*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - x[0]*x[0]*x[0]*x[1] + 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 7*x[0]*(x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]) + x[0]*x[0]*x[0]*(x[1]*x[1]) + 7*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 10*(x[1]*x[1]*x[1]) - 9*x[0]*(x[1]*x[1]*x[1]) - 4*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 7*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 3*x[2] + 8*x[0]*x[2] + x[0]*x[0]*x[2] + 2*(x[0]*x[0]*x[0])*x[2] - 3*(x[0]*x[0]*x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 3*x[1]*x[2] - 2*x[0]*x[1]*x[2] + 3*(x[0]*x[0])*x[1]*x[2] + 9*(x[0]*x[0]*x[0])*x[1]*x[2] + 5*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 7*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] - 3*(x[0]*x[0])*(x[1]*x[1])*x[2] + 6*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1]*x[1])*x[2] + x[0]*x[0]*(x[1]*x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1]*x[1])*x[2] + 10*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + x[2]*x[2] - 9*x[0]*(x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]) - 7*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) + 7*x[0]*x[1]*(x[2]*x[2]) - 2*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 5*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - x[1]*x[1]*(x[2]*x[2]) - 7*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 9*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - x[1]*x[1]*x[1]*(x[2]*x[2]) - 10*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - x[2]*x[2]*x[2] + 9*x[0]*(x[2]*x[2]*x[2]) - 10*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]) + 3*x[0]*x[1]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 3*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 7*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 8*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 4*(x[2]*x[2]*x[2]*x[2]) - 3*x[0]*(x[2]*x[2]*x[2]*x[2]) - 7*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + 6*x[1]*(x[2]*x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 4*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 9*(x[2]*x[2]*x[2]*x[2]*x[2]) - 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2])
        def dfunc(x):
            return [-7 - 4*x[0] - 9*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) - 50*(x[0]*x[0]*x[0]*x[0]) - 42*(x[0]*x[0]*x[0]*x[0]*x[0]) - 5*x[1] - 8*x[0]*x[1] - 3*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 45*(x[0]*x[0]*x[0]*x[0])*x[1] - 7*(x[1]*x[1]) + 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) + 28*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 9*(x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]) - 12*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]*x[1]) + 8*x[2] + 2*x[0]*x[2] + 6*(x[0]*x[0])*x[2] - 12*(x[0]*x[0]*x[0])*x[2] + 50*(x[0]*x[0]*x[0]*x[0])*x[2] - 2*x[1]*x[2] + 6*x[0]*x[1]*x[2] + 27*(x[0]*x[0])*x[1]*x[2] + 20*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 18*(x[0]*x[0])*(x[1]*x[1])*x[2] - 6*(x[1]*x[1]*x[1])*x[2] + 2*x[0]*(x[1]*x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1]*x[1])*x[2] - 9*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 21*(x[0]*x[0])*(x[2]*x[2]) + 32*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 7*x[1]*(x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]) + 15*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 7*(x[1]*x[1])*(x[2]*x[2]) - 18*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*(x[2]*x[2]*x[2]) - 20*x[0]*(x[2]*x[2]*x[2]) + 24*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 3*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) - 7*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 3*(x[2]*x[2]*x[2]*x[2]) - 14*x[0]*(x[2]*x[2]*x[2]*x[2]) - 4*x[1]*(x[2]*x[2]*x[2]*x[2]) - 10*(x[2]*x[2]*x[2]*x[2]*x[2]), 7 - 5*x[0] - 4*(x[0]*x[0]) - x[0]*x[0]*x[0] + 4*(x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 8*x[1] - 14*x[0]*x[1] + 18*(x[0]*x[0])*x[1] + 2*(x[0]*x[0]*x[0])*x[1] + 14*(x[0]*x[0]*x[0]*x[0])*x[1] - 30*(x[1]*x[1]) - 27*x[0]*(x[1]*x[1]) - 12*(x[0]*x[0])*(x[1]*x[1]) - 12*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 28*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 40*x[0]*(x[1]*x[1]*x[1]*x[1]) + 42*(x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] - 2*x[0]*x[2] + 3*(x[0]*x[0])*x[2] + 9*(x[0]*x[0]*x[0])*x[2] + 5*(x[0]*x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 6*(x[0]*x[0])*x[1]*x[2] + 12*(x[0]*x[0]*x[0])*x[1]*x[2] + 30*(x[1]*x[1])*x[2] - 18*x[0]*(x[1]*x[1])*x[2] + 3*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 40*x[0]*(x[1]*x[1]*x[1])*x[2] + 50*(x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) + 7*x[0]*(x[2]*x[2]) - 2*(x[0]*x[0])*(x[2]*x[2]) + 5*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 14*x[0]*x[1]*(x[2]*x[2]) - 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 3*(x[1]*x[1])*(x[2]*x[2]) - 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 12*(x[1]*x[1]*x[1])*(x[2]*x[2]) - x[2]*x[2]*x[2] + 3*x[0]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 6*x[1]*(x[2]*x[2]*x[2]) - 14*x[0]*x[1]*(x[2]*x[2]*x[2]) + 24*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]) - 4*x[0]*(x[2]*x[2]*x[2]*x[2]) - 8*x[1]*(x[2]*x[2]*x[2]*x[2]) - x[2]*x[2]*x[2]*x[2]*x[2], 3 + 8*x[0] + x[0]*x[0] + 2*(x[0]*x[0]*x[0]) - 3*(x[0]*x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] - 2*x[0]*x[1] + 3*(x[0]*x[0])*x[1] + 9*(x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0])*x[1] + 7*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 6*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 10*(x[1]*x[1]*x[1]) - 6*x[0]*(x[1]*x[1]*x[1]) + x[0]*x[0]*(x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]) + 10*x[0]*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*x[2] - 18*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 14*(x[0]*x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] + 14*x[0]*x[1]*x[2] - 4*(x[0]*x[0])*x[1]*x[2] + 10*(x[0]*x[0]*x[0])*x[1]*x[2] - 2*(x[1]*x[1])*x[2] - 14*x[0]*(x[1]*x[1])*x[2] - 18*(x[0]*x[0])*(x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1])*x[2] - 20*x[0]*(x[1]*x[1]*x[1])*x[2] + 6*(x[1]*x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 27*x[0]*(x[2]*x[2]) - 30*(x[0]*x[0])*(x[2]*x[2]) + 24*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 3*x[1]*(x[2]*x[2]) + 9*x[0]*x[1]*(x[2]*x[2]) - 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 9*(x[1]*x[1])*(x[2]*x[2]) - 21*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 16*(x[2]*x[2]*x[2]) - 12*x[0]*(x[2]*x[2]*x[2]) - 28*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 24*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) - 16*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 45*(x[2]*x[2]*x[2]*x[2]) - 50*x[0]*(x[2]*x[2]*x[2]*x[2]) - 5*x[1]*(x[2]*x[2]*x[2]*x[2]) - 30*(x[2]*x[2]*x[2]*x[2]*x[2])]
        def ddfunc(x):
            return [[-4 - 18*x[0] - 12*(x[0]*x[0]) - 200*(x[0]*x[0]*x[0]) - 210*(x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 6*x[0]*x[1] + 48*(x[0]*x[0])*x[1] - 180*(x[0]*x[0]*x[0])*x[1] + 18*(x[1]*x[1]) + 6*x[0]*(x[1]*x[1]) + 84*(x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 24*x[0]*(x[1]*x[1]*x[1]) - 18*(x[1]*x[1]*x[1]*x[1]) + 2*x[2] + 12*x[0]*x[2] - 36*(x[0]*x[0])*x[2] + 200*(x[0]*x[0]*x[0])*x[2] + 6*x[1]*x[2] + 54*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] + 36*x[0]*(x[1]*x[1])*x[2] + 2*(x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 42*x[0]*(x[2]*x[2]) + 96*(x[0]*x[0])*(x[2]*x[2]) - 4*x[1]*(x[2]*x[2]) + 30*x[0]*x[1]*(x[2]*x[2]) - 18*(x[1]*x[1])*(x[2]*x[2]) - 20*(x[2]*x[2]*x[2]) + 48*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 14*(x[2]*x[2]*x[2]*x[2]), -5 - 8*x[0] - 3*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) - 14*x[1] + 36*x[0]*x[1] + 6*(x[0]*x[0])*x[1] + 56*(x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) - 72*x[0]*(x[1]*x[1]*x[1]) - 40*(x[1]*x[1]*x[1]*x[1]) - 2*x[2] + 6*x[0]*x[2] + 27*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] + 8*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 36*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] + 6*x[0]*(x[1]*x[1])*x[2] + 40*(x[1]*x[1]*x[1])*x[2] + 7*(x[2]*x[2]) - 4*x[0]*(x[2]*x[2]) + 15*(x[0]*x[0])*(x[2]*x[2]) - 14*x[1]*(x[2]*x[2]) - 36*x[0]*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 14*x[1]*(x[2]*x[2]*x[2]) - 4*(x[2]*x[2]*x[2]*x[2]), 8 + 2*x[0] + 6*(x[0]*x[0]) - 12*(x[0]*x[0]*x[0]) + 50*(x[0]*x[0]*x[0]*x[0]) - 2*x[1] + 6*x[0]*x[1] + 27*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 6*x[0]*(x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) + 2*x[0]*(x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]) - 18*x[2] - 4*x[0]*x[2] - 42*(x[0]*x[0])*x[2] + 64*(x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 8*x[0]*x[1]*x[2] + 30*(x[0]*x[0])*x[1]*x[2] - 14*(x[1]*x[1])*x[2] - 36*x[0]*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 27*(x[2]*x[2]) - 60*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 21*(x[1]*x[1])*(x[2]*x[2]) - 12*(x[2]*x[2]*x[2]) - 56*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 50*(x[2]*x[2]*x[2]*x[2])], [-5 - 8*x[0] - 3*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) - 14*x[1] + 36*x[0]*x[1] + 6*(x[0]*x[0])*x[1] + 56*(x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) - 72*x[0]*(x[1]*x[1]*x[1]) - 40*(x[1]*x[1]*x[1]*x[1]) - 2*x[2] + 6*x[0]*x[2] + 27*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] + 8*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 36*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] + 6*x[0]*(x[1]*x[1])*x[2] + 40*(x[1]*x[1]*x[1])*x[2] + 7*(x[2]*x[2]) - 4*x[0]*(x[2]*x[2]) + 15*(x[0]*x[0])*(x[2]*x[2]) - 14*x[1]*(x[2]*x[2]) - 36*x[0]*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 14*x[1]*(x[2]*x[2]*x[2]) - 4*(x[2]*x[2]*x[2]*x[2]), 8 - 14*x[0] + 18*(x[0]*x[0]) + 2*(x[0]*x[0]*x[0]) + 14*(x[0]*x[0]*x[0]*x[0]) - 60*x[1] - 54*x[0]*x[1] - 24*(x[0]*x[0])*x[1] - 24*(x[0]*x[0]*x[0])*x[1] + 84*(x[1]*x[1]) - 48*x[0]*(x[1]*x[1]) - 108*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[1]*x[1]*x[1]) - 160*x[0]*(x[1]*x[1]*x[1]) + 210*(x[1]*x[1]*x[1]*x[1]) + 14*x[2] + 8*x[0]*x[2] - 6*(x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0])*x[2] + 60*x[1]*x[2] - 36*x[0]*x[1]*x[2] + 6*(x[0]*x[0])*x[1]*x[2] - 96*(x[1]*x[1])*x[2] + 120*x[0]*(x[1]*x[1])*x[2] + 200*(x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 14*x[0]*(x[2]*x[2]) - 18*(x[0]*x[0])*(x[2]*x[2]) - 6*x[1]*(x[2]*x[2]) - 60*x[0]*x[1]*(x[2]*x[2]) + 36*(x[1]*x[1])*(x[2]*x[2]) - 6*(x[2]*x[2]*x[2]) - 14*x[0]*(x[2]*x[2]*x[2]) + 48*x[1]*(x[2]*x[2]*x[2]) - 8*(x[2]*x[2]*x[2]*x[2]), -3 - 2*x[0] + 3*(x[0]*x[0]) + 9*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 14*x[1] + 8*x[0]*x[1] - 6*(x[0]*x[0])*x[1] + 12*(x[0]*x[0]*x[0])*x[1] + 30*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[1]*x[1]*x[1]) + 40*x[0]*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) - 4*x[2] + 14*x[0]*x[2] - 4*(x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 28*x[0]*x[1]*x[2] - 36*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 24*(x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 9*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) - 18*x[1]*(x[2]*x[2]) - 42*x[0]*x[1]*(x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 32*x[1]*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2])], [8 + 2*x[0] + 6*(x[0]*x[0]) - 12*(x[0]*x[0]*x[0]) + 50*(x[0]*x[0]*x[0]*x[0]) - 2*x[1] + 6*x[0]*x[1] + 27*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 6*x[0]*(x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) + 2*x[0]*(x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]) - 18*x[2] - 4*x[0]*x[2] - 42*(x[0]*x[0])*x[2] + 64*(x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 8*x[0]*x[1]*x[2] + 30*(x[0]*x[0])*x[1]*x[2] - 14*(x[1]*x[1])*x[2] - 36*x[0]*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 27*(x[2]*x[2]) - 60*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 21*(x[1]*x[1])*(x[2]*x[2]) - 12*(x[2]*x[2]*x[2]) - 56*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 50*(x[2]*x[2]*x[2]*x[2]), -3 - 2*x[0] + 3*(x[0]*x[0]) + 9*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 14*x[1] + 8*x[0]*x[1] - 6*(x[0]*x[0])*x[1] + 12*(x[0]*x[0]*x[0])*x[1] + 30*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[1]*x[1]*x[1]) + 40*x[0]*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) - 4*x[2] + 14*x[0]*x[2] - 4*(x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 28*x[0]*x[1]*x[2] - 36*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 24*(x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 9*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) - 18*x[1]*(x[2]*x[2]) - 42*x[0]*x[1]*(x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 32*x[1]*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]), 2 - 18*x[0] - 2*(x[0]*x[0]) - 14*(x[0]*x[0]*x[0]) + 16*(x[0]*x[0]*x[0]*x[0]) - 4*x[1] + 14*x[0]*x[1] - 4*(x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0])*x[1] - 2*(x[1]*x[1]) - 14*x[0]*(x[1]*x[1]) - 18*(x[0]*x[0])*(x[1]*x[1]) - 2*(x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]) - 6*x[2] + 54*x[0]*x[2] - 60*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] - 6*x[1]*x[2] + 18*x[0]*x[1]*x[2] - 48*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] - 42*x[0]*(x[1]*x[1])*x[2] + 48*(x[1]*x[1]*x[1])*x[2] + 48*(x[2]*x[2]) - 36*x[0]*(x[2]*x[2]) - 84*(x[0]*x[0])*(x[2]*x[2]) + 72*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 48*(x[1]*x[1])*(x[2]*x[2]) + 180*(x[2]*x[2]*x[2]) - 200*x[0]*(x[2]*x[2]*x[2]) - 20*x[1]*(x[2]*x[2]*x[2]) - 150*(x[2]*x[2]*x[2]*x[2])]]
elif funcType == "septic":
    if dimension == 1:
        def func(x):
            return 3 - 7*x[0] + 8*(x[0]*x[0]) + 7*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0])
        def dfunc(x):
            return [-7 + 16*x[0] + 21*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) - 40*(x[0]*x[0]*x[0]*x[0]) - 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 42*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])]
        def ddfunc(x):
            return [[16 + 42*x[0] - 60*(x[0]*x[0]) - 160*(x[0]*x[0]*x[0]) - 270*(x[0]*x[0]*x[0]*x[0]) - 252*(x[0]*x[0]*x[0]*x[0]*x[0])]]
    elif dimension == 2:
        def func(x):
            return -9 - 10*x[0] - 3*(x[0]*x[0]) - 10*(x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 6*x[1] + 8*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0])*x[1] - x[0]*x[0]*x[0]*x[0]*x[1] + 4*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 4*(x[1]*x[1]) + 10*x[0]*(x[1]*x[1]) - 8*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 9*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 8*(x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 6*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 6*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 8*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 6*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]*x[1]) + 9*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + 5*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 5*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1])
        def dfunc(x):
            return [-10 - 6*x[0] - 30*(x[0]*x[0]) - 24*(x[0]*x[0]*x[0]) - 30*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 56*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*x[1] - 8*x[0]*x[1] - 27*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0]*x[0])*x[1] - 24*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 10*(x[1]*x[1]) - 16*x[0]*(x[1]*x[1]) + 24*(x[0]*x[0])*(x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 40*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 4*(x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 24*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 12*x[0]*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]*x[1]) + 10*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 5*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]), -6 + 8*x[0] - 4*(x[0]*x[0]) - 9*(x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0] + 4*(x[0]*x[0]*x[0]*x[0]*x[0]) - 4*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 8*x[1] + 20*x[0]*x[1] - 16*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] + 18*(x[0]*x[0]*x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 24*(x[1]*x[1]) - 12*x[0]*(x[1]*x[1]) + 27*(x[0]*x[0])*(x[1]*x[1]) + 18*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 18*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 32*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 24*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 12*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 45*x[0]*(x[1]*x[1]*x[1]*x[1]) + 25*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 12*(x[1]*x[1]*x[1]*x[1]*x[1]) - 30*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 70*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1])]
        def ddfunc(x):
            return [[-6 - 60*x[0] - 72*(x[0]*x[0]) - 120*(x[0]*x[0]*x[0]) + 60*(x[0]*x[0]*x[0]*x[0]) + 336*(x[0]*x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 54*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 120*(x[0]*x[0]*x[0]*x[0])*x[1] - 16*(x[1]*x[1]) + 48*x[0]*(x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]) + 160*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 18*(x[1]*x[1]*x[1]) + 36*x[0]*(x[1]*x[1]*x[1]) + 72*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 12*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]), 8 - 8*x[0] - 27*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) - 24*(x[0]*x[0]*x[0]*x[0]*x[0]) + 20*x[1] - 32*x[0]*x[1] + 48*(x[0]*x[0])*x[1] + 72*(x[0]*x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0]*x[0])*x[1] - 12*(x[1]*x[1]) + 54*x[0]*(x[1]*x[1]) + 54*(x[0]*x[0])*(x[1]*x[1]) + 72*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) + 48*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 50*x[0]*(x[1]*x[1]*x[1]*x[1]) - 30*(x[1]*x[1]*x[1]*x[1]*x[1])], [8 - 8*x[0] - 27*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) - 24*(x[0]*x[0]*x[0]*x[0]*x[0]) + 20*x[1] - 32*x[0]*x[1] + 48*(x[0]*x[0])*x[1] + 72*(x[0]*x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0]*x[0])*x[1] - 12*(x[1]*x[1]) + 54*x[0]*(x[1]*x[1]) + 54*(x[0]*x[0])*(x[1]*x[1]) + 72*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) + 48*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 50*x[0]*(x[1]*x[1]*x[1]*x[1]) - 30*(x[1]*x[1]*x[1]*x[1]*x[1]), -8 + 20*x[0] - 16*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) + 18*(x[0]*x[0]*x[0]*x[0]) + 16*(x[0]*x[0]*x[0]*x[0]*x[0]) + 48*x[1] - 24*x[0]*x[1] + 54*(x[0]*x[0])*x[1] + 36*(x[0]*x[0]*x[0])*x[1] + 36*(x[0]*x[0]*x[0]*x[0])*x[1] + 96*(x[1]*x[1]) - 48*x[0]*(x[1]*x[1]) + 72*(x[0]*x[0])*(x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 180*(x[1]*x[1]*x[1]) + 180*x[0]*(x[1]*x[1]*x[1]) + 100*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 60*(x[1]*x[1]*x[1]*x[1]) - 150*x[0]*(x[1]*x[1]*x[1]*x[1]) - 420*(x[1]*x[1]*x[1]*x[1]*x[1])]]
    else:
        def func(x):
            return -10 - 4*x[0] + 10*(x[0]*x[0]) + 7*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0]*x[0] + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 9*x[1] + 9*x[0]*x[1] - x[0]*x[0]*x[1] - 7*(x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) - 10*x[0]*(x[1]*x[1]) - 10*(x[0]*x[0])*(x[1]*x[1]) - 7*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 5*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 10*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 9*x[0]*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 10*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1] + x[0]*(x[1]*x[1]*x[1]*x[1]) - 8*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]) + 6*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 2*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1]*x[1] - 6*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 3*x[2] - 3*x[0]*x[2] - 10*(x[0]*x[0])*x[2] - 6*(x[0]*x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0]*x[0])*x[2] + 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 3*x[1]*x[2] + 7*x[0]*x[1]*x[2] + 7*(x[0]*x[0])*x[1]*x[2] - 4*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + x[1]*x[1]*x[2] + 8*x[0]*(x[1]*x[1])*x[2] + 2*(x[0]*x[0])*(x[1]*x[1])*x[2] + 10*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 2*(x[1]*x[1]*x[1])*x[2] - 8*x[0]*(x[1]*x[1]*x[1])*x[2] - x[0]*x[0]*(x[1]*x[1]*x[1])*x[2] + 7*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 4*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + x[0]*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] - 9*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + 10*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 10*(x[0]*x[0])*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 10*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]) - x[0]*x[0]*x[1]*(x[2]*x[2]) - 2*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) + 8*(x[1]*x[1])*(x[2]*x[2]) + 2*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) + 4*(x[0]*x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 4*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 10*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 8*(x[0]*x[0])*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 4*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*x[0]*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 5*x[0]*(x[2]*x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]*x[2]) - 9*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - x[0]*x[0]*x[0]*x[0]*(x[2]*x[2]*x[2]) + 3*x[1]*(x[2]*x[2]*x[2]) + 2*x[0]*x[1]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 2*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 2*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 10*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - x[0]*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) - 6*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 7*(x[2]*x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) - 8*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + x[1]*(x[2]*x[2]*x[2]*x[2]) + x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - x[0]*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) + 8*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) - 9*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]*x[2]) + 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) + 7*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2])
        def dfunc(x):
            return [-4 + 20*x[0] + 21*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) + 48*(x[0]*x[0]*x[0]*x[0]*x[0]) + 56*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 9*x[1] - 2*x[0]*x[1] - 21*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 25*(x[0]*x[0]*x[0]*x[0])*x[1] + 54*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 10*(x[1]*x[1]) - 20*x[0]*(x[1]*x[1]) - 21*(x[0]*x[0])*(x[1]*x[1]) + 20*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 50*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 9*(x[1]*x[1]*x[1]) - 6*x[0]*(x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 40*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1] - 16*x[0]*(x[1]*x[1]*x[1]*x[1]) + 27*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 6*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] - 20*x[0]*x[2] - 18*(x[0]*x[0])*x[2] + 40*(x[0]*x[0]*x[0])*x[2] + 45*(x[0]*x[0]*x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 7*x[1]*x[2] + 14*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 8*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] + 30*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1])*x[2] - 2*x[0]*(x[1]*x[1]*x[1])*x[2] + 21*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 4*(x[1]*x[1]*x[1]*x[1])*x[2] + 18*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + x[1]*x[1]*x[1]*x[1]*x[1]*x[2] - 2*(x[2]*x[2]) - 20*x[0]*(x[2]*x[2]) - 12*(x[0]*x[0])*(x[2]*x[2]) + 24*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 50*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 8*x[1]*(x[2]*x[2]) - 2*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) + 2*(x[1]*x[1])*(x[2]*x[2]) - 20*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 12*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 16*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - 5*(x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]) - 27*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) + 2*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 10*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 2*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + x[1]*(x[2]*x[2]*x[2]*x[2]) - 2*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 9*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 10*(x[2]*x[2]*x[2]*x[2]*x[2]) + 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]), -9 + 9*x[0] - x[0]*x[0] - 7*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 16*x[1] - 20*x[0]*x[1] - 20*(x[0]*x[0])*x[1] - 14*(x[0]*x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0]*x[0])*x[1] - 20*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 24*(x[1]*x[1]) - 27*x[0]*(x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]) - 9*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 30*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) + 4*x[0]*(x[1]*x[1]*x[1]) - 32*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) + 30*x[0]*(x[1]*x[1]*x[1]*x[1]) - 10*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 36*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 56*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] + 7*x[0]*x[2] + 7*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] + 4*(x[0]*x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 2*x[1]*x[2] + 16*x[0]*x[1]*x[2] + 4*(x[0]*x[0])*x[1]*x[2] + 20*(x[0]*x[0]*x[0])*x[1]*x[2] - 16*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 6*(x[1]*x[1])*x[2] - 24*x[0]*(x[1]*x[1])*x[2] - 3*(x[0]*x[0])*(x[1]*x[1])*x[2] + 21*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 20*(x[1]*x[1]*x[1])*x[2] - 16*x[0]*(x[1]*x[1]*x[1])*x[2] + 36*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 10*(x[1]*x[1]*x[1]*x[1])*x[2] + 5*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] - 54*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 8*x[0]*(x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]) - 2*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 16*x[1]*(x[2]*x[2]) + 4*x[0]*x[1]*(x[2]*x[2]) - 20*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 12*(x[1]*x[1])*(x[2]*x[2]) - 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 16*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 36*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 50*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) + 2*x[0]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 2*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - 4*x[1]*(x[2]*x[2]*x[2]) + 20*x[0]*x[1]*(x[2]*x[2]*x[2]) - 2*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 18*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 24*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] + x[0]*(x[2]*x[2]*x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]*x[2]*x[2]) + 16*x[1]*(x[2]*x[2]*x[2]*x[2]) - 18*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) - x[2]*x[2]*x[2]*x[2]*x[2] + 5*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 2*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 7*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]), 3 - 3*x[0] - 10*(x[0]*x[0]) - 6*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] + 7*x[0]*x[1] + 7*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] + 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + x[1]*x[1] + 8*x[0]*(x[1]*x[1]) + 2*(x[0]*x[0])*(x[1]*x[1]) + 10*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 2*(x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]) - x[0]*x[0]*(x[1]*x[1]*x[1]) + 7*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 5*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) - 2*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 9*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 20*x[2] - 4*x[0]*x[2] - 20*(x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 16*x[0]*x[1]*x[2] - 2*(x[0]*x[0])*x[1]*x[2] - 4*(x[0]*x[0]*x[0])*x[1]*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 16*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] - 20*(x[0]*x[0])*(x[1]*x[1])*x[2] + 8*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1])*x[2] - 20*x[0]*(x[1]*x[1]*x[1])*x[2] + 16*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1]*x[1])*x[2] + 18*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) - 15*x[0]*(x[2]*x[2]) - 3*(x[0]*x[0])*(x[2]*x[2]) - 27*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 3*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) + 6*x[0]*x[1]*(x[2]*x[2]) - 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 6*(x[1]*x[1])*(x[2]*x[2]) + 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 3*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 18*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 18*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 28*(x[2]*x[2]*x[2]) - 8*x[0]*(x[2]*x[2]*x[2]) - 32*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 32*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) + 4*x[0]*x[1]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 32*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 36*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 24*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 30*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) - 5*x[1]*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 5*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 36*(x[2]*x[2]*x[2]*x[2]*x[2]) - 12*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 42*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 14*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2])]
        def ddfunc(x):
            return [[20 + 42*x[0] - 60*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) + 240*(x[0]*x[0]*x[0]*x[0]) + 336*(x[0]*x[0]*x[0]*x[0]*x[0]) - 2*x[1] - 42*x[0]*x[1] + 60*(x[0]*x[0])*x[1] + 100*(x[0]*x[0]*x[0])*x[1] + 270*(x[0]*x[0]*x[0]*x[0])*x[1] - 20*(x[1]*x[1]) - 42*x[0]*(x[1]*x[1]) + 60*(x[0]*x[0])*(x[1]*x[1]) - 200*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]) + 120*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 16*(x[1]*x[1]*x[1]*x[1]) + 54*x[0]*(x[1]*x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]*x[1]) - 20*x[2] - 36*x[0]*x[2] + 120*(x[0]*x[0])*x[2] + 180*(x[0]*x[0]*x[0])*x[2] + 60*(x[0]*x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 24*x[0]*x[1]*x[2] + 48*(x[0]*x[0])*x[1]*x[2] - 160*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] + 60*x[0]*(x[1]*x[1])*x[2] - 96*(x[0]*x[0])*(x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1])*x[2] + 42*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 20*(x[2]*x[2]) - 24*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 200*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 12*x[0]*x[1]*(x[2]*x[2]) - 48*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 20*(x[1]*x[1])*(x[2]*x[2]) + 24*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 16*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 2*(x[2]*x[2]*x[2]) - 54*x[0]*(x[2]*x[2]*x[2]) - 12*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) + 12*x[0]*x[1]*(x[2]*x[2]*x[2]) - 2*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 16*(x[2]*x[2]*x[2]*x[2]) - 48*x[0]*(x[2]*x[2]*x[2]*x[2]) - 2*x[1]*(x[2]*x[2]*x[2]*x[2]) + 10*(x[2]*x[2]*x[2]*x[2]*x[2]), 9 - 2*x[0] - 21*(x[0]*x[0]) + 20*(x[0]*x[0]*x[0]) + 25*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 20*x[1] - 40*x[0]*x[1] - 42*(x[0]*x[0])*x[1] + 40*(x[0]*x[0]*x[0])*x[1] - 100*(x[0]*x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) - 64*x[0]*(x[1]*x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) - 36*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*x[2] + 14*x[0]*x[2] - 12*(x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0])*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[2] + 16*x[1]*x[2] + 8*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 64*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 63*(x[0]*x[0])*(x[1]*x[1])*x[2] - 16*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 8*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 6*(x[0]*x[0])*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 4*x[1]*(x[2]*x[2]) - 40*x[0]*x[1]*(x[2]*x[2]) + 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 48*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 36*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 2*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 20*x[1]*(x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 18*x[1]*(x[2]*x[2]*x[2]*x[2]) + 5*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 - 20*x[0] - 18*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 45*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] + 14*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 40*(x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) + 30*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 2*x[0]*(x[1]*x[1]*x[1]) + 21*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1] - 4*x[2] - 40*x[0]*x[2] - 24*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] + 100*(x[0]*x[0]*x[0]*x[0])*x[2] - 16*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] - 32*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 40*x[0]*(x[1]*x[1])*x[2] + 24*(x[0]*x[0])*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 32*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 15*(x[2]*x[2]) - 6*x[0]*(x[2]*x[2]) - 81*(x[0]*x[0])*(x[2]*x[2]) - 12*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) + 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 30*(x[1]*x[1])*(x[2]*x[2]) - 6*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 8*(x[2]*x[2]*x[2]) - 64*x[0]*(x[2]*x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]*x[2]) - 36*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 50*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*x[1]*(x[2]*x[2]*x[2]*x[2]) - 12*(x[2]*x[2]*x[2]*x[2]*x[2])], [9 - 2*x[0] - 21*(x[0]*x[0]) + 20*(x[0]*x[0]*x[0]) + 25*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 20*x[1] - 40*x[0]*x[1] - 42*(x[0]*x[0])*x[1] + 40*(x[0]*x[0]*x[0])*x[1] - 100*(x[0]*x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) - 64*x[0]*(x[1]*x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) - 36*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*x[2] + 14*x[0]*x[2] - 12*(x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0])*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[2] + 16*x[1]*x[2] + 8*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 64*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 63*(x[0]*x[0])*(x[1]*x[1])*x[2] - 16*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 8*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 6*(x[0]*x[0])*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 4*x[1]*(x[2]*x[2]) - 40*x[0]*x[1]*(x[2]*x[2]) + 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 48*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 36*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 2*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 20*x[1]*(x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 18*x[1]*(x[2]*x[2]*x[2]*x[2]) + 5*(x[2]*x[2]*x[2]*x[2]*x[2]), 16 - 20*x[0] - 20*(x[0]*x[0]) - 14*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) - 20*(x[0]*x[0]*x[0]*x[0]*x[0]) - 48*x[1] - 54*x[0]*x[1] - 18*(x[0]*x[0])*x[1] - 18*(x[0]*x[0]*x[0])*x[1] + 60*(x[0]*x[0]*x[0]*x[0])*x[1] + 12*(x[1]*x[1]) + 12*x[0]*(x[1]*x[1]) - 96*(x[0]*x[0])*(x[1]*x[1]) + 108*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 200*(x[1]*x[1]*x[1]) + 120*x[0]*(x[1]*x[1]*x[1]) - 40*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 180*x[0]*(x[1]*x[1]*x[1]*x[1]) - 336*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*x[2] + 16*x[0]*x[2] + 4*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] - 16*(x[0]*x[0]*x[0]*x[0])*x[2] + 12*x[1]*x[2] - 48*x[0]*x[1]*x[2] - 6*(x[0]*x[0])*x[1]*x[2] + 42*(x[0]*x[0]*x[0])*x[1]*x[2] + 60*(x[1]*x[1])*x[2] - 48*x[0]*(x[1]*x[1])*x[2] + 108*(x[0]*x[0])*(x[1]*x[1])*x[2] - 40*(x[1]*x[1]*x[1])*x[2] + 20*x[0]*(x[1]*x[1]*x[1])*x[2] - 270*(x[1]*x[1]*x[1]*x[1])*x[2] + 16*(x[2]*x[2]) + 4*x[0]*(x[2]*x[2]) - 20*(x[0]*x[0])*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 24*x[1]*(x[2]*x[2]) - 60*x[0]*x[1]*(x[2]*x[2]) + 48*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 48*(x[1]*x[1])*(x[2]*x[2]) + 108*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 200*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 4*(x[2]*x[2]*x[2]) + 20*x[0]*(x[2]*x[2]*x[2]) - 2*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 36*x[1]*(x[2]*x[2]*x[2]) + 36*x[0]*x[1]*(x[2]*x[2]*x[2]) - 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 16*(x[2]*x[2]*x[2]*x[2]) - 18*x[0]*(x[2]*x[2]*x[2]*x[2]) + 36*x[1]*(x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 + 7*x[0] + 7*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 4*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*x[1] + 16*x[0]*x[1] + 4*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 21*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 20*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) + 5*x[0]*(x[1]*x[1]*x[1]*x[1]) - 54*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[2] - 16*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[2] + 32*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 40*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 48*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] - 100*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) + 6*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 12*x[1]*(x[2]*x[2]) + 60*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 54*(x[1]*x[1])*(x[2]*x[2]) + 54*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 4*(x[2]*x[2]*x[2]) + 4*x[0]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 64*x[1]*(x[2]*x[2]*x[2]) - 72*x[0]*x[1]*(x[2]*x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*(x[2]*x[2]*x[2]*x[2]) - 10*x[1]*(x[2]*x[2]*x[2]*x[2]) + 42*(x[2]*x[2]*x[2]*x[2]*x[2])], [-3 - 20*x[0] - 18*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 45*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] + 14*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 40*(x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) + 30*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 2*x[0]*(x[1]*x[1]*x[1]) + 21*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1] - 4*x[2] - 40*x[0]*x[2] - 24*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] + 100*(x[0]*x[0]*x[0]*x[0])*x[2] - 16*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] - 32*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 40*x[0]*(x[1]*x[1])*x[2] + 24*(x[0]*x[0])*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 32*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 15*(x[2]*x[2]) - 6*x[0]*(x[2]*x[2]) - 81*(x[0]*x[0])*(x[2]*x[2]) - 12*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) + 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 30*(x[1]*x[1])*(x[2]*x[2]) - 6*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 8*(x[2]*x[2]*x[2]) - 64*x[0]*(x[2]*x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]*x[2]) - 36*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 50*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*x[1]*(x[2]*x[2]*x[2]*x[2]) - 12*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 + 7*x[0] + 7*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 4*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*x[1] + 16*x[0]*x[1] + 4*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 21*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 20*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) + 5*x[0]*(x[1]*x[1]*x[1]*x[1]) - 54*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[2] - 16*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[2] + 32*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 40*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 48*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] - 100*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) + 6*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 12*x[1]*(x[2]*x[2]) + 60*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 54*(x[1]*x[1])*(x[2]*x[2]) + 54*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 4*(x[2]*x[2]*x[2]) + 4*x[0]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 64*x[1]*(x[2]*x[2]*x[2]) - 72*x[0]*x[1]*(x[2]*x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*(x[2]*x[2]*x[2]*x[2]) - 10*x[1]*(x[2]*x[2]*x[2]*x[2]) + 42*(x[2]*x[2]*x[2]*x[2]*x[2]), 20 - 4*x[0] - 20*(x[0]*x[0]) - 8*(x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]*x[0]) - 4*x[1] - 16*x[0]*x[1] - 2*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] - 8*(x[0]*x[0]*x[0]*x[0])*x[1] + 16*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) - 20*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]) + 16*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1]*x[1]) + 18*x[2] - 30*x[0]*x[2] - 6*(x[0]*x[0])*x[2] - 54*(x[0]*x[0]*x[0])*x[2] - 6*(x[0]*x[0]*x[0]*x[0])*x[2] + 18*x[1]*x[2] + 12*x[0]*x[1]*x[2] - 48*(x[0]*x[0])*x[1]*x[2] + 12*(x[0]*x[0]*x[0])*x[1]*x[2] - 12*(x[1]*x[1])*x[2] + 60*x[0]*(x[1]*x[1])*x[2] - 6*(x[0]*x[0])*(x[1]*x[1])*x[2] - 36*(x[1]*x[1]*x[1])*x[2] + 36*x[0]*(x[1]*x[1]*x[1])*x[2] - 36*(x[1]*x[1]*x[1]*x[1])*x[2] + 84*(x[2]*x[2]) - 24*x[0]*(x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]) - 96*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 12*x[1]*(x[2]*x[2]) + 12*x[0]*x[1]*(x[2]*x[2]) - 12*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 96*(x[1]*x[1])*(x[2]*x[2]) - 108*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 120*(x[2]*x[2]*x[2]) + 200*x[0]*(x[2]*x[2]*x[2]) + 100*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 20*x[1]*(x[2]*x[2]*x[2]) + 100*x[0]*x[1]*(x[2]*x[2]*x[2]) - 20*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 180*(x[2]*x[2]*x[2]*x[2]) - 60*x[0]*(x[2]*x[2]*x[2]*x[2]) + 210*x[1]*(x[2]*x[2]*x[2]*x[2]) - 84*(x[2]*x[2]*x[2]*x[2]*x[2])]]
elif funcType == "sinusoidal":
    if dimension == 1:
        def func(x):
            return 3. + 2.*np.sin(5*x[0])
        def dfunc(x):
            return [10.*np.cos(5*x[0])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0])]]
    elif dimension == 2:
        def func(x):
            return 3. + 2.*(np.sin(5*x[0]) + np.sin(5*x[1]))
        def dfunc(x):
            return [10.*np.cos(5*x[0]), 10.*np.cos(5*x[1])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0]), 0], [0, -50.*np.sin(5*x[1])]]
    else:
        def func(x):
            return 3. + 2.*(np.sin(5*x[0]) + np.sin(5*x[1]) + np.sin(5*x[2]))
        def dfunc(x):
            return [10.*np.cos(5*x[0]), 10.*np.cos(5*x[1]), 10.*np.cos(5*x[2])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0]), 0, 0], [0, -50.*np.sin(5*x[1]), 0], [0, 0, -50.*np.sin(5*x[2])]]
else:
    raise ValueError, "function type {} not found".format(funcType)

#-------------------------------------------------------------------------------
# Create RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(dataBase = dataBase,
                           W = WT,
                           volumeType = volumeType,
                           needHessian = testHessian)
packages = [rk]

#-------------------------------------------------------------------------------
# Create a state directly and initialize physics package
#-------------------------------------------------------------------------------
connectivity = dataBase.connectivityMap()
state = State(dataBase, packages)
derivs = StateDerivatives(dataBase, packages)
rk.initializeProblemStartup(dataBase)
rk.registerState(dataBase, state)
rk.registerDerivatives(dataBase, derivs)
rk.preStepInitialize(dataBase, state, derivs)

#-------------------------------------------------------------------------------
# Get data from state
#-------------------------------------------------------------------------------
position = state.vectorFields(HydroFieldNames.position)
H = state.symTensorFields(HydroFieldNames.H)
volume = state.scalarFields(HydroFieldNames.volume)
zerothCorrections = state.vector_of_doubleFields(HydroFieldNames.rkZerothCorrections)
corrections = state.vector_of_doubleFields(HydroFieldNames.rkCorrections)
# normal = state.vectorFields(HydroFieldNames.normal)

#-------------------------------------------------------------------------------
# Compute corrections
#-------------------------------------------------------------------------------
if computeCorrectionsDirectly:
    rk_time = time.time()
    RKUtilities.computeCorrections(connectivity, WT, volume, position, H, testHessian,
                                   zerothCorrections, corrections)
    rk_time = time.time() - rk_time
else:
    rk_time = time.time()
    rk.initialize(0.0, 0.0, dataBase, state, derivs)
    rk_time = time.time() - rk_time
output("rk_time")

#-------------------------------------------------------------------------------
# Get old corrections
#-------------------------------------------------------------------------------
if correctionOrder <= QuadraticOrder and checkAgainstOld:
    A = dataBase.newFluidScalarFieldList(name="A")
    B = dataBase.newFluidVectorFieldList(name="B")
    C = dataBase.newFluidTensorFieldList(name="C")
    gradA = dataBase.newFluidVectorFieldList(name="gradA")
    gradB = dataBase.newFluidTensorFieldList(name="gradB")
    gradC = dataBase.newFluidThirdRankTensorFieldList(name="gradB")
    
    M0 = dataBase.newFluidScalarFieldList(name="M0")
    M1 = dataBase.newFluidVectorFieldList(name="M1")
    M2 = dataBase.newFluidSymTensorFieldList(name="M2")
    M3 = dataBase.newFluidThirdRankTensorFieldList(name="M3")
    M4 = dataBase.newFluidFourthRankTensorFieldList(name="M4")
    gradM0 = dataBase.newFluidVectorFieldList(name="grad M0")
    gradM1 = dataBase.newFluidTensorFieldList(name="grad M1")
    gradM2 = dataBase.newFluidThirdRankTensorFieldList(name="grad M2")
    gradM3 = dataBase.newFluidFourthRankTensorFieldList(name="grad M3")
    gradM4 = dataBase.newFluidFifthRankTensorFieldList(name="grad M4")
    
    surfacePoint = dataBase.newFluidIntFieldList(name="surface point")

    old_rk_time = time.time()
    computeCRKSPHMoments(connectivity, WT, volume, position, H, correctionOrder, NodeCoupling(),
                         M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4)
    computeCRKSPHCorrections(M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4,
                             H, surfacePoint, correctionOrder,
                             A, B, C, gradA, gradB, gradC)
    old_rk_time = time.time() - old_rk_time
    output("old_rk_time")

if quitAfterTiming:
    quit()
    
#-------------------------------------------------------------------------------
# Set up a simple method to calculate the kernel
#-------------------------------------------------------------------------------
# Base kernel
def getBaseKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hj = H(nj, j)
    w = RKUtilities.evaluateBaseKernel(WT, xij, Hj)
    dw = RKUtilities.evaluateBaseGradient(WT, xij, Hj)
    if testHessian:
        ddw = RKUtilities.evaluateBaseHessian(WT, xij, Hj)
        ddw = np.reshape(ddw, (dimension, dimension))
    else:
        ddw = np.zeros((dimension, dimension))
    return w, dw, ddw

# RK kernel
def getNewKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hj = H(nj, j)
    c = corrections(ni, i)
    w = RKUtilities.evaluateKernel(WT, xij, Hj, c)
    dw = RKUtilities.evaluateGradient(WT, xij, Hj, c)
    if testHessian:
        ddwtemp = RKUtilities.evaluateHessian(WT, xij, Hj, c)
        ddw = np.zeros((dimension, dimension))
        for d1 in range(dimension):
            for d2 in range(dimension):
                ddw[d1,d2] = ddwtemp(d1, d2)
    else:
        ddw = np.zeros((dimension, dimension))
    return w, dw, ddw

# Old kernel
def getOldKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    Hj = H(nj, j)
    xij = xi - xj
    etaj = Hj * xij
    Hdetj = Hj.Determinant()
    Ai = A(ni, i)
    Bi = B(ni, i) if correctionOrder >= LinearOrder else Vector.zero
    Ci = C(ni, i) if correctionOrder >= QuadraticOrder else Tensor.zero
    dAi = gradA(ni, i)
    dBi = gradB(ni, i) if correctionOrder >= LinearOrder else Tensor.zero
    dCi = gradC(ni, i) if correctionOrder >= QuadraticOrder else ThirdRankTensor.zero
    w, dwtemp, dw = CRKSPHKernelAndGradient(WT, correctionOrder,
                                            xij, etaj, Hj, Hdetj,
                                            Ai, Bi, Ci, dAi, dBi, dCi)
    ddw = Tensor.zero
    return w, dw, ddw

# Kernel to test with
getKernel = getOldKernel if useOldKernel else getNewKernel
getKernel = getBaseKernel if useBaseKernel else getNewKernel

#-------------------------------------------------------------------------------
# Check against old method of doing corrections
#-------------------------------------------------------------------------------
if correctionOrder <= QuadraticOrder and checkAgainstOld:
    # Check old and new corrections
    for i in range(nodes.numNodes):
        ni = 0
        # Get old corrections in same format as new
        old = [A(ni, i)]
        if correctionOrder >= LinearOrder:
            for d in range(dimension):
                old.append(A(ni,i) * B(ni, i)[d])
        if correctionOrder >= QuadraticOrder:
            for d1 in range(dimension):
                for d2 in range(d1, dimension):
                    old.append(A(ni, i) * C(ni, i)[dimension*d1+d2])
        for d in range(dimension):
            old.append(gradA(ni, i)[d])
        if correctionOrder >= LinearOrder:
            a = A(ni, i)
            for d2 in range(dimension):
                da = gradA(ni, i)[d2]
                for d1 in range(dimension):
                    b = B(ni, i)[d1]
                    db = gradB(ni, i)[d1 * dimension + d2]
                    old.append(da * b + a * db)
        if correctionOrder >= QuadraticOrder:
            for d1 in range(dimension):
                for d2 in range(d1, dimension):
                    c = C(ni, i)[d1*dimension + d2]
                    for d3 in range(dimension):
                        da = gradA(ni, i)[d3]
                        dc = gradC(ni, i)[d1*dimension*dimension + d2*dimension + d3]
                        old.append(da * c + a * dc)
        # Compare to new corrections
        c = np.array(corrections(ni, i))
        err = np.abs(np.subtract(c, old))
        if any(err > 1.e-6) and printErrors:
            print i
            print "\t", c
            print "\t", old
            print "\t", np.subtract(c, old)

    # Compare kernel values
    for i in range(nodes.numNodes):
        ni = 0
        connectivityi = connectivity.connectivityForNode(ni, i)
        for nj, neighbors in enumerate(connectivityi):
            for j in neighbors:
                wnew, dwnew, ddwnew = getNewKernel(ni, i, nj, j)
                wold, dwold, ddwtemp = getOldKernel(ni, i, nj, j)
                dwnew = np.array(dwnew)
                dwold = np.array(dwold)
                werr = np.abs(wnew - wold)
                dwerr = np.abs(dwnew - dwold)
                if werr > 1.e-6 or any(dwerr > 1.e-6) and printErrors:
                    print i, j, wnew, wold, dwnew, dwold

#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
# Get the nodes to check
if numToCheck == -1:
    nodesToCheck = range(nodes.numNodes)
elif numToCheck > 0:
    nodesToCheck = random.sample(range(nodes.numNodes), numToCheck)
else:
    raise ValueError, "numToCheck must be -1 or positive"

interp_time = time.time()
vals = np.zeros((nodes.numNodes, 2))
dvals = np.zeros((nodes.numNodes, dimension, 2))
ddvals = np.zeros((nodes.numNodes, dimension, dimension, 2))
ni = 0
for i in nodesToCheck:#range(nodes.numNodes):
    xi = position(ni, i)
    fi = func(xi)
    def addToValues(nj, j):
        xj = position(nj, j)
        if type(xj) is not type(xi):
            raise TypeError, "error in xj, i = {}, j = {}".format(i, j)
        fj = func(xj)
        xij = xi - xj
        vj = volume(nj, j)
        w, dw, ddw = getKernel(ni, i, nj, j)
        vals[i,0] += vj * w * fj
        dvals[i,:,0] += vj * dw * fj
        if testHessian:
            ddvals[i,:,:,0] += vj * ddw * fj
    connectivityi = connectivity.connectivityForNode(ni, i)
    for nj, neighbors in enumerate(connectivityi):
        for j in neighbors:
            addToValues(nj, j)
    addToValues(ni, i)
    vals[i,1] = fi
    dvals[i,:,1] = dfunc(xi)
    if testHessian:
        ddvals[i,:,:,1] = ddfunc(xi)

interp_time = time.time() - interp_time
output("interp_time")

#-------------------------------------------------------------------------------
# Optionally plot results
#-------------------------------------------------------------------------------
def plotThings(numLocal, anaLocal, title):
    num = mpi.gather(numLocal)
    ana = mpi.gather(anaLocal)
    if mpi.rank == 0:
        num = np.concatenate(num)
        ana = np.concatenate(ana)
        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
        fig, ax1 = plt.subplots()
        ax1.plot(num, label="num", color=colors[0], linestyle="-")
        ax1.plot(ana, label="ana", color=colors[1], linestyle="--")
        plt.legend()
        ax2 = ax1.twinx()
        err = np.divide(num - ana, np.mean(np.abs(ana)) + 1.e-15)
        ax2.semilogy(err, label="err", color=colors[2], linestyle="-.")
        plt.title(title)
        plt.legend()
if plot:
    plotThings(vals[:,0], vals[:,1], "vals")
    for d in range(dimension):
        plotThings(dvals[:,d,0], dvals[:,d,1], "dvals{}".format(d))
    if testHessian:
        for d1 in range(dimension):
            for d2 in range(dimension):
                plotThings(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1], "ddvals{}-{}".format(d1, d2))
    plt.show()
    
#-------------------------------------------------------------------------------
# Check error
#-------------------------------------------------------------------------------
def getError(num, ana):
    err = np.sum(np.abs(num - ana))
    tot = np.sum(np.abs(ana))
    err = mpi.allreduce(err, mpi.SUM)
    tot = mpi.allreduce(tot, mpi.SUM)
    return err / (tot + 1.e-11)
error = getError(vals[:,0], vals[:,1])
output("error")
derror = [getError(dvals[:,d,0], dvals[:,d,1]) for d in range(dimension)]
output("derror")
if testHessian:
    dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
    output("dderror")

if correctionOrder == RKOrder.ZerothOrder:
    ni = 0
    zerothErr = np.zeros((nodes.numNodes))
    for i in range(nodes.numNodes):
        zerothErr[i] = getError(np.array(corrections(ni, i)), np.array(zerothCorrections(ni, i)))
    output("np.amax(zerothErr)")
    if any([e0 > tolerance for e0 in zerothErr]):
        raise ValueError, "zeroth corrections do not agree"

if error > tolerance:
    raise ValueError, "error is greater than tolerance"
if funcType != "constant":
    if any([de > 10*tolerance for de in derror]):
        raise ValueError, "gradient error is greater than tolerance"
if testHessian and funcType != "constant" and funcType != "linear":
    if any([dde > 100*tolerance for dde in dderror]):
        raise ValueError, "hessian error is greater than tolerance"
        
