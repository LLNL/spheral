#ATS:test(SELF, "--dimension 1 --correctionOrder ZerothOrder --funcType constant", label="RK zeroth 1d")
#ATS:test(SELF, "--dimension 1 --correctionOrder LinearOrder --funcType linear", label="RK linear 1d")
#ATS:test(SELF, "--dimension 1 --correctionOrder QuadraticOrder --funcType quadratic", label="RK quadratic 1d")
#ATS:test(SELF, "--dimension 1 --correctionOrder CubicOrder --funcType cubic", label="RK cubic 1d")
#ATS:test(SELF, "--dimension 2 --correctionOrder ZerothOrder --funcType constant", label="RK zeroth 2d")
#ATS:test(SELF, "--dimension 2 --correctionOrder LinearOrder --funcType linear", label="RK linear 2d")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuadraticOrder --funcType quadratic", label="RK quadratic 2d")
#ATS:test(SELF, "--dimension 2 --correctionOrder CubicOrder --funcType cubic", label="RK cubic 2d")
#ATS:test(SELF, "--dimension 3 --correctionOrder ZerothOrder --funcType constant", label="RK zeroth 3d")
#ATS:test(SELF, "--dimension 3 --correctionOrder LinearOrder --funcType linear", label="RK linear 3d")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuadraticOrder --funcType quadratic", label="RK quadratic 3d")
#ATS:test(SELF, "--dimension 3 --correctionOrder CubicOrder --funcType cubic", label="RK cubic 3d")

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
    correctionOrder = LinearOrder,
    needHessian = True,
    
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
elif dimension == 2:
    from Spheral2d import *
else:
    from Spheral3d import *

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

#-------------------------------------------------------------------------------
# Build the RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(dataBase = dataBase,
                   W = WT,
                   correctionOrder = correctionOrder,
                   volumeType = CRKMassOverDensity,
                   needHessian = needHessian)
packages = [rk]

#-------------------------------------------------------------------------------
# Create integrator
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(dataBase)
for p in packages:
    integrator.appendPhysicsPackage(p)
for p in packages:
    p.initializeProblemStartup(dataBase)
        
#-------------------------------------------------------------------------------
# Get a state and apply the applicable methods directly to it
#-------------------------------------------------------------------------------
dataBase.updateConnectivityMap(True)

state = State(dataBase, integrator.physicsPackages())
derivs = StateDerivatives(dataBase, integrator.physicsPackages())
integrator.preStepInitialize(state, derivs)

rk_time = time.time()
integrator.initializeDerivatives(0.0, 0.0, state, derivs)
rk_time = time.time() - rk_time

integrator.applyGhostBoundaries(state, derivs)
integrator.finalizeGhostBoundaries()

output("rk_time")

#-------------------------------------------------------------------------------
# Make sure changes to H propagate to corrections
#-------------------------------------------------------------------------------
# if initagain:
#     dataBase.updateConnectivityMap(True)
#     rk.initializeProblemStartup(dataBase)


# integrator.applyGhostBoundaries(state, derivs)
# integrator.finalizeGhostBoundaries()

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
# Get some data
#-------------------------------------------------------------------------------
position = state.vectorFields(HydroFieldNames.position)
H = state.symTensorFields(HydroFieldNames.H)
volume = state.scalarFields(HydroFieldNames.volume)

A = state.scalarFields(HydroFieldNames.A_RK)
dA = state.vectorFields(HydroFieldNames.gradA_RK)
ddA = state.tensorFields(HydroFieldNames.hessA_RK)
B = state.vectorFields(HydroFieldNames.B_RK)
dB = state.tensorFields(HydroFieldNames.gradB_RK)
ddB = state.thirdRankTensorFields(HydroFieldNames.hessB_RK)
C = state.tensorFields(HydroFieldNames.C_RK)
dC = state.thirdRankTensorFields(HydroFieldNames.gradC_RK)
ddC = state.fourthRankTensorFields(HydroFieldNames.hessC_RK)
D = state.thirdRankTensorFields(HydroFieldNames.D_RK)
dD = state.fourthRankTensorFields(HydroFieldNames.gradD_RK)
ddD = state.fifthRankTensorFields(HydroFieldNames.hessD_RK)

connectivity = dataBase.connectivityMap()

#-------------------------------------------------------------------------------
# Get zeroth-order correction to check against
#-------------------------------------------------------------------------------
# A_check = np.zeros(nodes.numNodes)
# if correctionOrder == ZerothOrder:
#     ni = 0
#     nj = 0
#     for i in range(nodes.numNodes):
#         xi = position(ni, i)
#         m0 = 0.
#         connectivityi = np.append(connectivity.connectivityForNode(ni, i), i)
#         for j in connectivityi:
#             xj = position(nj, j)
#             xij = xi - xj
#             Hij = H(nj, j)
#             etaij = Hij * xij
#             wij = WT(etaij, Hij)
#             vj = volume(nj, j)
#             m0 += vj * wij
#         A_check[i] = 1. / m0

#-------------------------------------------------------------------------------
# Set up a simple method to calculate the kernel
#-------------------------------------------------------------------------------
interp_time = time.time()
b = Vector.zero
db = Tensor.zero
ddb = ThirdRankTensor.zero

dxij = Tensor.zero
dxij.Identity()

zeroTensor = Tensor.zero

# Base kernel
def getBaseKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hij = H(nj, j)
    etaij = Hij * xij
    H2ij = Hij.square()
    etaijmag = etaij.magnitude()
    etaijmaginv = 0.0 if etaijmag < 1.e-13 else 1. / etaijmag 
    Hetaij = Hij * etaij * etaijmaginv
    Heta2ij = Hetaij.selfdyad()
    w = WT(etaij, Hij)
    dw = Hij * etaij * etaijmaginv * WT.grad(etaij, Hij)
    ddw = Tensor((H2ij - Heta2ij) * etaijmaginv * WT.grad(etaij, Hij) + Heta2ij * WT.grad2(etaij, Hij))
    return w, dw, ddw

# RK kernel
def getKernel(ni, i, nj, j):
    xi = position(ni, i)
    a = A(ni, i)
    da = dA(ni, i)
    dda = ddA(ni, i)
    if correctionOrder >=LinearOrder:
        b = B(ni, i)
        db = dB(ni, i)
        ddb = ddB(ni, i)
    else:
        b = Vector.zero
        db = Tensor.zero
        ddb = ThirdRankTensor.zero
    if correctionOrder >= QuadraticOrder:
        c = C(ni, i)
        dc = dC(ni, i)
        ddc = ddC(ni, i)
    else:
        c = Tensor.zero
        dc = ThirdRankTensor.zero
        ddc = FourthRankTensor.zero
    if correctionOrder >= CubicOrder:
        d = D(ni, i)
        dd = dD(ni, i)
        ddd = ddD(ni, i)
    else:
        d = ThirdRankTensor.zero
        dd = FourthRankTensor.zero
        ddd = FifthRankTensor.zero
    vj = volume(nj, j)
    xj = position(nj, j)
    xij = xi - xj
    Hij = H(nj, j)
    etaij = Hij * xij
    w = evaluateRKKernel(WT, correctionOrder,
                         etaij, Hij, xij, 
                         a, b, c, d)
    dw = evaluateRKGradient(WT, correctionOrder,
                            etaij, Hij, xij, dxij,
                            a, b, c, d,
                            da, db, dc, dd)
    if needHessian:
        ddw = evaluateRKHessian(WT, correctionOrder,
                                etaij, Hij, xij, dxij,
                                a, b, c, d,
                                da, db, dc, dd,
                                dda, ddb, ddc, ddd)
        ddw = np.reshape(ddw, (dimension, dimension))
    else:
        ddw = zeroTensor
    return w, dw, ddw

#-------------------------------------------------------------------------------
# Check 2d linear against analytic solutions
#-------------------------------------------------------------------------------
if dimension == 2 and correctionOrder == LinearOrder and checkConditions:
    ni = 0
    m0 = np.zeros((nodes.numNodes))
    m1 = np.zeros((nodes.numNodes, dimension))
    m2 = np.zeros((nodes.numNodes, dimension, dimension))
    dm0 = np.zeros((nodes.numNodes, dimension))
    dm1 = np.zeros((nodes.numNodes, dimension, dimension))
    dm2 = np.zeros((nodes.numNodes, dimension, dimension, dimension))
    ddm0 = np.zeros((nodes.numNodes, dimension, dimension))
    ddm1 = np.zeros((nodes.numNodes, dimension, dimension, dimension))
    ddm2 = np.zeros((nodes.numNodes, dimension, dimension, dimension, dimension))
    aPy = np.zeros((nodes.numNodes))
    bPy = np.zeros((nodes.numNodes, dimension))
    daPy = np.zeros((nodes.numNodes, dimension))
    dbPy = np.zeros((nodes.numNodes, dimension, dimension))
    ddaPy = np.zeros((nodes.numNodes, dimension, dimension))
    ddbPy = np.zeros((nodes.numNodes, dimension, dimension, dimension))

    # Compute m values
    for i in range(nodes.numNodes):
        connectivityi = connectivity.connectivityForNode(ni, i)
        xi = position(ni, i)
        def addToValues(nj, j):
            xj = position(nj, j)
            xij = xi - xj
            dxij = np.identity(dimension)
            vj = volume(nj, j)
            w, dw, ddw = getBaseKernel(ni, i, nj, j)
            m0[i] += vj * w
            m1[i] += vj * xij * w
            m2[i] += vj * np.outer(xij, xij) * w
            for k1 in range(dimension):
                dm0[i,k1] += vj * dw[k1]
                dm1[i,k1] += vj * (dxij[k1,:] * w + xij * dw[k1])
                dm2[i,k1] += vj * (np.outer(xij, xij) * dw[k1] + (np.outer(dxij[:,k1], xij) + np.outer(xij, dxij[:,k1])) * w)
            return
        addToValues(ni,i)
        for nj, neighbors in enumerate(connectivityi):
            for j in neighbors:
                addToValues(nj, j)
    # Compute corrections
    for i in range(nodes.numNodes):
        mmat = np.array([[m0[i],   m1[i,0],   m1[i,1]],
                         [m1[i,0], m2[i,0,0], m2[i,0,1]],
                         [m1[i,1], m2[i,1,0], m2[i,1,1]]])
        rhs = np.array([1., 0., 0.])
        coeff = np.linalg.solve(mmat, rhs)
        aPy[i] = coeff[0]
        bPy[i,0] = coeff[1]
        bPy[i,1] = coeff[2]
        for k1 in range(dimension):
            dmmat = np.array([[dm0[i,k1],   dm1[i,k1,0],   dm1[i,k1,1]],
                              [dm1[i,k1,0], dm2[i,k1,0,0], dm2[i,k1,0,1]],
                              [dm1[i,k1,1], dm2[i,k1,1,0], dm2[i,k1,1,1]]])
            rhs = -np.matmul(dmmat, coeff)
            dcoeff = np.linalg.solve(mmat, rhs)
            daPy[i,k1] = dcoeff[0]
            dbPy[i,k1,0] = dcoeff[1]
            dbPy[i,k1,1] = dcoeff[2]
            
    # Check corrections
    ni = 0
    for i in range(nodes.numNodes):
        diffA = A(ni, i) - aPy[i]
        diffB = B(ni, i) - bPy[i]
        diffdA = dA(ni, i) - daPy[i]
        diffdB = np.reshape(dB(ni, i), (2,2)) - dbPy[i]

    def getKernelPy(ni, i, nj, j):
        xi = position(ni, i)
        a = aPy[i]
        da = Vector.zero
        dda = Tensor.zero
        b = Vector.zero
        db = Tensor.zero
        ddb = ThirdRankTensor.zero
        for k1 in range(dimension):
            da[k1] = daPy[i,k1]
            b[k1] = bPy[i,k1]
            for k2 in range(dimension):
                dda[2*k1+k2] = ddaPy[i,k2,k1]
                db[2*k1+k2] = dbPy[i,k2,k1]
                for k3 in range(dimension):
                    ddb[4*k1+2*k2+k3] = ddbPy[i,k3,k2,k1]
        c = Tensor.zero
        dc = ThirdRankTensor.zero
        ddc = FourthRankTensor.zero
        d = ThirdRankTensor.zero
        dd = FourthRankTensor.zero
        ddd = FifthRankTensor.zero
        vj = volume(nj, j)
        xj = position(nj, j)
        xij = xi - xj
        Hij = H(nj, j)
        etaij = Hij * xij
        w, dw, ddw = getBaseKernel(ni, i, nj, j)
        wr = (a + np.dot(xij, b)) * w
        dwr = (da + b + np.multiply(xij[0], [db[0], db[1]]) + np.multiply(xij[1], [db[2],db[3]])) * w + (a + np.dot(xij, b)) * dw
        # w = evaluateRKKernel(WT, correctionOrder,
        #                      etaij, Hij, xij, 
        #                      a, b, c, d)
        # dw = evaluateRKGradient(WT, correctionOrder,
        #                         etaij, Hij, xij, dxij,
        #                         a, b, c, d,
        #                         da, db, dc, dd)
        # if needHessian:
        #     ddw = evaluateRKHessian(WT, correctionOrder,
        #                             etaij, Hij, xij, dxij,
        #                             a, b, c, d,
        #                             da, db, dc, dd,
        #                             dda, ddb, ddc, ddd)
        #     ddw = np.reshape(ddw, (dimension, dimension))
        # else:
        ddwr = zeroTensor
        return wr, dwr, ddwr
    # getKernel = getKernelPy
# import code
# code.interact(local=locals())
#-------------------------------------------------------------------------------
# Make sure linear conditions are satisfied
#-------------------------------------------------------------------------------
if correctionOrder >= LinearOrder and checkConditions:
    val0 = np.zeros((nodes.numNodes))
    val1 = np.zeros((nodes.numNodes, dimension))
    dval0 = np.zeros((nodes.numNodes, dimension))
    dval1 = np.zeros((nodes.numNodes, dimension, dimension))
    ddval0 = np.zeros((nodes.numNodes, dimension, dimension))
    ddval1 = np.zeros((nodes.numNodes, dimension, dimension, dimension))
    ni = 0
    for i in range(nodes.numNodes):
        connectivityi = connectivity.connectivityForNode(ni, i)
        xi = position(ni, i)
        def addToValues(nj, j):
            xj = position(nj, j)
            xij = xi - xj
            dxij = np.identity(dimension)
            vj = volume(nj, j)
            w, dw, ddw = getKernel(ni, i, nj, j)
            val0[i] += vj * w
            val1[i,:] += vj * xij * w
            dval0[i,:] += vj * dw
            dval1[i,:] += vj * (np.outer(xij, dw) + dxij * w)
            return
        addToValues(ni, i)
        for nj, neighbors in enumerate(connectivityi):
            for j in neighbors:
                addToValues(nj, j)
    val0max = np.amax(np.abs(val0 - 1))
    val1max = np.amax(np.abs(val1))
    dval0max = np.amax(np.abs(dval0))
    dval1max = np.amax(np.abs(dval1))
    output("val0max,val1max,dval0max,dval1max")
# quit()
            
#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
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
        if needHessian:
            ddvals[i,:,:,0] += vj * ddw * fj
    connectivityi = connectivity.connectivityForNode(ni, i)
    for nj, neighbors in enumerate(connectivityi):
        for j in neighbors:
            addToValues(nj, j)
    addToValues(ni, i)
    vals[i,1] = fi
    dvals[i,:,1] = dfunc(xi)
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
        ax2.plot(err, label="err", color=colors[2], linestyle="-.")
        plt.title(title)
        plt.legend()
if plot:
    plotThings(vals[:,0], vals[:,1], "vals")
    for d in range(dimension):
        plotThings(dvals[:,d,0], dvals[:,d,1], "dvals{}".format(d))
    if needHessian:
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
derror = [getError(dvals[:,d,0], dvals[:,d,1]) for d in range(dimension)]
dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
derrormax = np.amax(derror)
output("error")
output("derror")
output("dderror")
if error > tolerance:
    raise ValueError, "error is greater than tolerance"
