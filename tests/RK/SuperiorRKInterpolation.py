#-------------------------------------------------------------------------------
# Manufactured diffusion test
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

    # CRK options
    correctionOrder = LinearOrder,
    volumeType = CRKMassOverDensity,
    testHessian = False,
    useOldKernel = False, # Test using old kernel
    
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
    tolerance = 1.e-14,
    
    # Plotting
    plot = False,
)

if dimension == 1:
    from Spheral1d import *
    if correctionOrder == ZerothOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d0
        SuperiorRKUtilities = SuperiorRKUtilities1d0
    elif correctionOrder == LinearOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d1
        SuperiorRKUtilities = SuperiorRKUtilities1d1
    elif correctionOrder == QuadraticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d2
        SuperiorRKUtilities = SuperiorRKUtilities1d2
    elif correctionOrder == CubicOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d3
        SuperiorRKUtilities = SuperiorRKUtilities1d3
    elif correctionOrder == QuarticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d4
        SuperiorRKUtilities = SuperiorRKUtilities1d4
    elif correctionOrder == QuinticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections1d5
        SuperiorRKUtilities = SuperiorRKUtilities1d5
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)
elif dimension == 2:
    from Spheral2d import *
    if correctionOrder == ZerothOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d0
        SuperiorRKUtilities = SuperiorRKUtilities2d0
    elif correctionOrder == LinearOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d1
        SuperiorRKUtilities = SuperiorRKUtilities2d1
    elif correctionOrder == QuadraticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d2
        SuperiorRKUtilities = SuperiorRKUtilities2d2
    elif correctionOrder == CubicOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d3
        SuperiorRKUtilities = SuperiorRKUtilities2d3
    elif correctionOrder == QuarticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d4
        SuperiorRKUtilities = SuperiorRKUtilities2d4
    elif correctionOrder == QuinticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections2d5
        SuperiorRKUtilities = SuperiorRKUtilities2d5
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)
else:
    from Spheral3d import *
    if correctionOrder == ZerothOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d0
        SuperiorRKUtilities = SuperiorRKUtilities3d0
    elif correctionOrder == LinearOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d1
        SuperiorRKUtilities = SuperiorRKUtilities3d1
    elif correctionOrder == QuadraticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d2
        SuperiorRKUtilities = SuperiorRKUtilities3d2
    elif correctionOrder == CubicOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d3
        SuperiorRKUtilities = SuperiorRKUtilities3d3
    elif correctionOrder == QuarticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d4
        SuperiorRKUtilities = SuperiorRKUtilities3d4
    elif correctionOrder == QuinticOrder:
        SuperiorRKCorrections = SuperiorRKCorrections3d5
        SuperiorRKUtilities = SuperiorRKUtilities3d5
    else:
        raise ValueError, "correction order \"{}\" not found".format(correctionOrder)

# if mpi.procs > 1:
#     raise ValueError, "need to add parallel boundaries and error calculation"
    
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
    def func(x):
        return 2.0 + 3.0 * np.sum(x)
    def dfunc(x):
        return 3.0 * np.ones_like(x)
    def ddfunc(x):
        return [[0.0 for v in x] for v in x]
elif funcType == "quadratic":
    if dimension == 1:
        def func(x):
            return 2 + 3*x[0] + 4*np.power(x[0],2)
        def dfunc(x):
            return [3 + 8*x[0]]
        def ddfunc(x):
            return [[8]]
    elif dimension == 2:
        def func(x):
            return 2 + 3*(x[0] + x[1]) + 4*(np.power(x[0],2) + x[0]*x[1] + np.power(x[1],2))
        def dfunc(x):
            return [3 + 4*(2*x[0] + x[1]), 3 + 4*(x[0] + 2*x[1])]
        def ddfunc(x):
            return [[8, 4], [4, 8]]
    else:
        def func(x):
            return 2 + 3*(x[0] + x[1] + x[2]) + 4*(np.power(x[0],2) + x[0]*x[1] + np.power(x[1],2) + x[0]*x[2] + x[1]*x[2] + np.power(x[2],2))
        def dfunc(x):
            return [3 + 4*(2*x[0] + x[1] + x[2]), 3 + 4*(x[0] + 2*x[1] + x[2]), 3 + 4*(x[0] + x[1] + 2*x[2])]
        def ddfunc(x):
            return [[8, 4, 4], [4, 8, 4], [4, 4, 8]]
elif funcType == "cubic":
    if dimension == 1:
        def func(x):
            return 2 + 3*x[0] + 4*np.power(x[0],2) + 5*np.power(x[0],3)
        def dfunc(x):
            return [3 + 8*x[0] + 15*np.power(x[0],2)]
        def ddfunc(x):
            return [[8 + 30*x[0]]]
    elif dimension == 2:
        def func(x):
            return 2 + 3*(x[0] + x[1]) + 4*(np.power(x[0],2) + x[0]*x[1] + np.power(x[1],2)) + 5*(np.power(x[0],3) + np.power(x[0],2)*x[1] + x[0]*np.power(x[1],2) + np.power(x[1],3))
        def dfunc(x):
            return [3 + 4*(2*x[0] + x[1]) + 5*(3*np.power(x[0],2) + 2*x[0]*x[1] + np.power(x[1],2)), 3 + 4*(x[0] + 2*x[1]) + 5*(np.power(x[0],2) + 2*x[0]*x[1] + 3*np.power(x[1],2))]
        def ddfunc(x):
            return [[8 + 5*(6*x[0] + 2*x[1]), 4 + 5*(2*x[0] + 2*x[1])], [4 + 5*(2*x[0] + 2*x[1]), 8 + 5*(2*x[0] + 6*x[1])]]
    else:
        def func(x):
            return 2 + 3*(x[0] + x[1] + x[2]) + 4*(np.power(x[0],2) + x[0]*x[1] + np.power(x[1],2) + x[0]*x[2] + x[1]*x[2] + np.power(x[2],2)) + 5*(np.power(x[0],3) + np.power(x[0],2)*x[1] + x[0]*np.power(x[1],2) + np.power(x[1],3) + np.power(x[0],2)*x[2] + x[0]*x[1]*x[2] + np.power(x[1],2)*x[2] + x[0]*np.power(x[2],2) + x[1]*np.power(x[2],2) + np.power(x[2],3))
        def dfunc(x):
            return [3 + 4*(2*x[0] + x[1] + x[2]) + 5*(3*np.power(x[0],2) + 2*x[0]*x[1] + np.power(x[1],2) + 2*x[0]*x[2] + x[1]*x[2] + np.power(x[2],2)), 3 + 4*(x[0] + 2*x[1] + x[2]) + 5*(np.power(x[0],2) + 2*x[0]*x[1] + 3*np.power(x[1],2) + x[0]*x[2] + 2*x[1]*x[2] + np.power(x[2],2)), 3 + 4*(x[0] + x[1] + 2*x[2]) + 5*(np.power(x[0],2) + x[0]*x[1] + np.power(x[1],2) + 2*x[0]*x[2] + 2*x[1]*x[2] + 3*np.power(x[2],2))]
        def ddfunc(x):
            return [[8 + 5*(6*x[0] + 2*x[1] + 2*x[2]), 4 + 5*(2*x[0] + 2*x[1] + x[2]), 4 + 5*(2*x[0] + x[1] + 2*x[2])], [4 + 5*(2*x[0] + 2*x[1] + x[2]), 8 + 5*(2*x[0] + 6*x[1] + 2*x[2]), 4 + 5*(x[0] + 2*x[1] + 2*x[2])], [4 + 5*(2*x[0] + x[1] + 2*x[2]), 4 + 5*(x[0] + 2*x[1] + 2*x[2]), 8 + 5*(2*x[0] + 2*x[1] + 6*x[2])]]
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
rk = SuperiorRKCorrections(dataBase = dataBase,
                           W = WT,
                           volumeType = volumeType,
                           needHessian = testHessian)
packages = [rk]

#-------------------------------------------------------------------------------
# Create a state directly and apply the applicable methods to it
#-------------------------------------------------------------------------------
connectivity = dataBase.connectivityMap()
state = State(dataBase, packages)
derivs = StateDerivatives(dataBase, packages)

rk.initializeProblemStartup(dataBase)
rk.registerState(dataBase, state)
rk.registerDerivatives(dataBase, derivs)
rk.preStepInitialize(dataBase, state, derivs)
rk_time = time.time()
rk.initialize(0.0, 0.0, dataBase, state, derivs)
rk_time = time.time() - rk_time
output("rk_time")

#-------------------------------------------------------------------------------
# Get data
#-------------------------------------------------------------------------------
position = state.vectorFields(HydroFieldNames.position)
H = state.symTensorFields(HydroFieldNames.H)
volume = state.scalarFields(HydroFieldNames.volume)
corrections = state.vector_of_doubleFields(HydroFieldNames.rkCorrections)

# for i in range(nodes.numNodes):
#     ni = 0
#     print i, corrections(ni, i)

#-------------------------------------------------------------------------------
# Get old corrections
#-------------------------------------------------------------------------------
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

computeCRKSPHMoments(connectivity, WT, volume, position, H, correctionOrder, NodeCoupling(),
                     M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4)
computeCRKSPHCorrections(M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4,
                         H, correctionOrder,
                         A, B, C, gradA, gradB, gradC)

#-------------------------------------------------------------------------------
# Set up a simple method to calculate the kernel
#-------------------------------------------------------------------------------
# Base kernel
def getBaseKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hj = H(nj, j)
    w = SuperiorRKUtilities.evaluateBaseKernel(WT, xij, Hj)
    dw = SuperiorRKUtilities.evaluateBaseGradient(WT, xij, Hj)
    if testHessian:
        ddw = SuperiorRKUtilities.evaluateBaseHessian(WT, xij, Hj)
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
    w = SuperiorRKUtilities.evaluateKernel(WT, xij, Hj, c)
    dw = SuperiorRKUtilities.evaluateGradient(WT, xij, Hj, c)
    if testHessian:
        ddw = SuperiorRKUtilities.evaluateHessian(WT, xij, Hj, c)
        ddw = np.reshape(ddw, (dimension, dimension))
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

#-------------------------------------------------------------------------------
# Check against old method of doing corrections
#-------------------------------------------------------------------------------
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
    if any(err > 1.e-6):
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
            if werr > 1.e-6 or any(dwerr > 1.e-6):
                print i, j, wnew, wold, dwnew, dwold

#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
interp_time = time.time()
vals = np.zeros((nodes.numNodes, 2))
dvals = np.zeros((nodes.numNodes, dimension, 2))
ddvals = np.zeros((nodes.numNodes, dimension, dimension, 2))
ni = 0
for i in range(nodes.numNodes):
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
    return err / (tot + 1.e-15)
error = getError(vals[:,0], vals[:,1])
output("error")
derror = [getError(dvals[:,d,0], dvals[:,d,1]) for d in range(dimension)]
output("derror")
if testHessian:
    dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
    output("dderror")
    
if error > tolerance:
    raise ValueError, "error is greater than tolerance"
if any([de > tolerance for de in derror]):
    raise ValueError, "gradient error is greater than tolerance"
if testHessian:
    if any([dde > tolerance for dde in dderror]):
        raise ValueError, "hessian error is greater than tolerance"
        
