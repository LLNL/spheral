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

if mpi.procs > 1:
    raise ValueError, "need to add parallel boundaries and error calculation"
    
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

dataBase.updateConnectivityMap(True)
connectivity = dataBase.connectivityMap()

#-------------------------------------------------------------------------------
# Build the RK object
#-------------------------------------------------------------------------------
rk_time = time.time()
rk = RKCorrections(dataBase = dataBase,
                   W = WT,
                   correctionOrder = correctionOrder,
                   volumeType = CRKMassOverDensity,
                   needHessian = needHessian)
packages = [rk]

#-------------------------------------------------------------------------------
# Run the startup stuff 
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(dataBase)
for p in packages:
    integrator.appendPhysicsPackage(p)
control = SpheralController(integrator, WT)

#-------------------------------------------------------------------------------
# Make sure changes to H propagate to corrections
#-------------------------------------------------------------------------------
rk.initializeProblemStartup(dataBase)
rk_time = time.time() - rk_time
output("rk_time")


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
position = dataBase.fluidPosition
H = dataBase.fluidHfield
volume = rk.volume
        
A = rk.A
dA = rk.gradA
ddA = rk.hessA
B = rk.B
dB = rk.gradB
ddB = rk.hessB
C = rk.C
dC = rk.gradC
ddC = rk.hessC
D = rk.D
dD = rk.gradD
ddD = rk.hessD

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
# Try interpolation
#-------------------------------------------------------------------------------
interp_time = time.time()
b = Vector.zero
db = Tensor.zero
ddb = ThirdRankTensor.zero
c = Tensor.zero
dc = ThirdRankTensor.zero
ddc = FourthRankTensor.zero
d = ThirdRankTensor.zero
dd = FourthRankTensor.zero
ddd = FifthRankTensor.zero

dxij = Tensor.zero
dxij.Identity()

vals = np.zeros((nodes.numNodes, 2))
dvals = np.zeros((nodes.numNodes, dimension, 2))
ddvals = np.zeros((nodes.numNodes, dimension, dimension, 2))
ni = 0
nj = 0

for i in range(nodes.numNodes):
    xi = position(ni, i)
    a = A(ni, i)
    da = dA(ni, i)
    dda = ddA(ni, i)
    if correctionOrder >=LinearOrder:
        b = B(ni, i)
        db = dB(ni, i)
        ddb = ddB(ni, i)
        if correctionOrder >= QuadraticOrder:
            c = C(ni, i)
            dc = dC(ni, i)
            ddc = ddC(ni, i)
            if correctionOrder >= CubicOrder:
                d = D(ni, i)
                dd = dD(ni, i)
                ddd = ddD(ni, i)
    connectivityi = np.append(connectivity.connectivityForNode(ni, i), i)
    print(c)
    fi = func(xi)
    for j in connectivityi:
        xj = position(nj, j)
        fj = func(xj)
        xij = xi - xj
        Hij = H(nj, j)
        etaij = Hij * xij
        vj = volume(nj, j)
        w = evaluateRKKernel(WT, correctionOrder,
                             etaij, Hij, xij, 
                             a, b, c, d)
        dw = evaluateRKGradient(WT, correctionOrder,
                                etaij, Hij, xij, dxij,
                                a, b, c, d,
                                da, db, dc, dd)
        vals[i,0] += vj * w * fj
        dvals[i,:,0] += vj * dw * (fj - fi)
        if needHessian:
            ddw = evaluateRKHessian(WT, correctionOrder,
                                    etaij, Hij, xij, dxij,
                                    a, b, c, d,
                                    da, db, dc, dd,
                                    dda, ddb, ddc, ddd)
            ddw = np.reshape(ddw, (dimension, dimension))
            ddvals[i,:,:,0] += vj * ddw * (fj - fi)
    vals[i,1] = fi
    dvals[i,:,1] = dfunc(xi)
    ddvals[i,:,:,1] = ddfunc(xi)

interp_time = time.time() - interp_time
output("interp_time")

#-------------------------------------------------------------------------------
# Optionally plot results
#-------------------------------------------------------------------------------
def plotThings(num, ana, title):
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
    return np.sum(np.abs(num - ana)) / (np.sum(np.abs(ana)) + 1.e-15)
error = getError(vals[:,0], vals[:,1])
derror = [getError(dvals[:,d,0], dvals[:,d,1]) for d in range(dimension)]
dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
derrormax = np.amax(derror)
output("error")
output("derror")
if error > tolerance:
    raise ValueError, "error is greater than tolerance"
