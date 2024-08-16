#ATS:t1 = test(SELF, "--dimension 1 --order 100 --tolerance 5.0e-3", label="integration, 1d", np=1)
#ATS:t2 = test(SELF, "--dimension 2 --nx 10 --ny 10 --order 10 --tolerance 4.0e-4", label="integration, 2d", np=1)
#ATS:t3 = test(SELF, "--dimension 3 --nx 5 --ny 5 --nz 5 --order 6", label="integration, 3d", np=1)
#ATS:r1 = test(SELF, "--dimension 1 --nx 20 --order 100 --correctionOrderIntegration 1", label="integration, 1d, rk1", np=1)
#ATS:r1 = test(SELF, "--dimension 1 --nx 20 --nPerh 10.01 --order 100 --correctionOrderIntegration 4", label="integration, 1d, rk4", np=1)
#ATS:r2 = test(SELF, "--dimension 2 --nx 20 --ny 20 --order 10 --correctionOrderIntegration 1", label="integration, 2d, rk1", np=1)
#ATS:t3 = test(SELF, "--dimension 3 --nx 5 --ny 5 --nz 5 --order 6 --correctionOrderIntegration 1", label="integration, 3d, rk1", np=1)

#-------------------------------------------------------------------------------
# Test the kernel integration methods
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
import time
from matplotlib import pyplot as plt
title("Kernel integration test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial parameters
    dimension = 1,
    nx = 20,
    ny = 20,
    nz = 20,
    x0 = -2.0,
    x1 = 2.0,
    y0 = -2.0,
    y1 = 2.0,
    z0 = -2.0,
    z1 = 2.0,

    # RK options
    useRK = False,
    correctionOrderIntegration = -1,
    volumeType = RKVoronoiVolume,

    # Integration parameters
    clipBoundaries = True,
    
    # Testing options
    randomizeNodes = False,
    ranfrac = 0.2,
    testHessian = False,
    printErrors = False,
    outputAllSurfaceVolume = False,
    checkAllSurfaceVolume = False,
    
    # Manufactured parameters
    funcType = "linear",
    
    # Integrator parameters
    order = 4,
    hMult = 1.0,
    
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
    dataDirBase = "dumps-KernelIntegration",

    # Error
    tolerance = 1.e-4,
    
    # Plotting
    plot = False,
)
correctionOrder = LinearOrder # We don't actually use this
useOverlap = False

if nPerh < correctionOrderIntegration:
    print("nPerh is not large enough for correction order: {} < {}".format(nPerh, correctionOrderIntegration))
    
if mpi.procs > 1:
    raise ValueError("test is not written for parallel")
    
#-------------------------------------------------------------------------------
# Choose correct dimension aliases
#-------------------------------------------------------------------------------
exec("from Spheral%id import *" % dimension)

#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------
# Set limits
lims = [[x0, x1], [y0, y1], [z0, z1]]
delta = [(x1 - x0)/nx, (y1 - y0) / ny, (z1 - z0) / nz]
length = [x1 - x0, y1 - y0, z1 - z0]
midpoints = [0.5*(x0+x1), 0.5*(y0+y1), 0.5*(z0+z1)]
analytic_volume = 1.
for d in range(dimension):
    analytic_volume *= lims[d][1] - lims[d][0]
if dimension == 1:
    analytic_surface = 2.
elif dimension == 2:
    analytic_surface = 2 * (length[0] + length[1])
elif dimension == 3:
    analytic_surface = 2 * (length[0] * length[1] +
                            length[0] * length[2] +
                            length[1] * length[2])
# Get point spacing
units = MKS()
deltaMin = min(delta[:dimension])
deltaMax = max(delta[:dimension])
hmin = deltaMin * nPerh * hminmult
hmax = deltaMax * nPerh * hmaxmult

output("delta")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "correctionOrder={}".format(correctionOrderIntegration),
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
WT = TableKernel(WendlandC2Kernel(), 1000)
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
    dx = delta[0]
    dy = delta[1]
    dz = delta[2]
    pos = nodes.positions()
    for i in range(nodes.numInternalNodes):
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
method = SPHSmoothingScale(IdealH, WT)
iterateIdealH(dataBase,
              [method],
              [],
              100, # max h iterations
              1.e-4) # h tolerance
dataBase.updateConnectivityMap(True, useOverlap) # need ghost and overlap connectivity

#-------------------------------------------------------------------------------
# Create RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(orders = set([ZerothOrder, correctionOrder]), 
                   dataBase = dataBase,
                   W = WT,
                   volumeType = volumeType,
                   needHessian = testHessian,
                   updateInFinalize = False)
output("rk")
output("rk.correctionOrders")
WR = rk.WR(correctionOrder)
output("WR")
packages = [rk]

#-------------------------------------------------------------------------------
# Add faceted boundaries
#-------------------------------------------------------------------------------
if clipBoundaries:
    if dimension == 1:
        points = vector_of_Vector([Vector(x0), Vector(x1)])
        facetedBoundary = Box1d(points)
    elif dimension == 2:
        points = vector_of_Vector([Vector(x0, y0), Vector(x0, y1),
                                   Vector(x1, y0), Vector(x1, y1)])
        facetedBoundary = Polygon(points)
    else:
        points = vector_of_Vector([Vector(x0, y0, z0), Vector(x0, y0, z1),
                                   Vector(x0, y1, z0), Vector(x0, y1, z1),
                                   Vector(x1, y0, z0), Vector(x1, y0, z1),
                                   Vector(x1, y1, z0), Vector(x1, y1, z1)])
        facetedBoundary = Polyhedron(points)
    output("facetedBoundary")
    output("facetedBoundary.volume")
    rk.addFacetedBoundary(facetedBoundary)

#-------------------------------------------------------------------------------
# Create a state directly and initialize physics package
#-------------------------------------------------------------------------------
connectivity = dataBase.connectivityMap(True, useOverlap)
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
zerothCorrections = state.RKCoefficientsFields(RKFieldNames.rkCorrections(ZerothOrder))
corrections = state.RKCoefficientsFields(RKFieldNames.rkCorrections(correctionOrder))
cells = state.facetedVolumeFields(HydroFieldNames.cells)
flags = state.vector_of_CellFaceFlagFields(HydroFieldNames.cellFaceFlags)
# normal = state.vectorFields(HydroFieldNames.normal)

#-------------------------------------------------------------------------------
# Compute corrections
#-------------------------------------------------------------------------------
WR.computeCorrections(connectivity, volume, position, H, testHessian,
                      zerothCorrections, corrections)

#-------------------------------------------------------------------------------
# Create flat connectivity
#-------------------------------------------------------------------------------
connectivity_time = time.time()
flatConnectivity = FlatConnectivity()
flatConnectivity.computeIndices(dataBase)
if useOverlap:
    flatConnectivity.computeOverlapIndices(dataBase)
flatConnectivity.computeSurfaceIndices(dataBase, state)
connectivity_time = time.time() - connectivity_time
output("connectivity_time")

#-------------------------------------------------------------------------------
# Create integrator, add integrals, and perform integration
#-------------------------------------------------------------------------------
if correctionOrderIntegration < 0:
    integrationKernel = SPHIntegrationKernel(WT)
else:
    integrationKernel = eval("RKIntegrationKernel{}d{}(WT)".format(dimension, correctionOrderIntegration))
output("integrationKernel")
integrator = KernelIntegrator(order, integrationKernel, dataBase, flatConnectivity)
vlK_f = LinearKernel()
vlG_f = LinearGrad()
vbKK_f = BilinearKernelKernel()
vbGK_f = BilinearGradKernel()
vbKG_f = BilinearKernelGrad()
vbGdG_f = BilinearGradDotGrad()
vbGpG_f = BilinearGradProdGrad()
slKn_f = LinearSurfaceNormalKernel()
sbKKn_f = BilinearSurfaceNormalKernelKernel()
sbKGdn_f = BilinearSurfaceNormalKernelDotGrad()
sbKKn2_f = BilinearSurfaceNormalKernelKernelFromGrad()
vcc_f = CellCoefficient()
scn_f = SurfaceNormalCoefficient()
volumeIntegrals = [vlK_f, vlG_f, vbKK_f, vbGK_f, vbKG_f, vbGdG_f, vbGpG_f, sbKKn2_f, vcc_f]
surfaceIntegrals = [slKn_f, sbKKn_f, sbKGdn_f, scn_f]
integrals = volumeIntegrals + surfaceIntegrals

for integral in integrals:
    integrator.addIntegral(integral)

integrator.setState(0.0, # time
                    state)

integration_time = time.time()
integrator.performIntegration()
integration_time = time.time() - integration_time
output("integration_time")
output("integrator.totalNumSubcells()")
output("integrator.totalNumSubfacets()")

# Get the integrals
vlK = list(vlK_f.values())
vlG = list(vlG_f.values())
vbKK = list(vbKK_f.values())
vbGK = list(vbGK_f.values())
vbKG = list(vbKG_f.values())
vbGdG = list(vbGdG_f.values())
vbGpG = list(vbGpG_f.values())
slKn = list(slKn_f.values())
sbKKn = list(sbKKn_f.values())
sbKGdn = list(sbKGdn_f.values())
sbKKn2 = list(sbKKn2_f.values())
slKn2 = list(vlG_f.values()) # The surface integral is the linear integral for constant coefficients
vcc = list(vcc_f.values())
scn = list(scn_f.values())

#-------------------------------------------------------------------------------
# Check volumes
#-------------------------------------------------------------------------------
checksum = 0

totcell = 0.
totvol = 0.
totint = 0.
for i in range(nodes.numNodes):
    totcell += cells[0][i].volume
    totvol += volume[0][i]
    totint += vcc[i]
    # print cells[0][i].xmin.x, cells[0][i].xmax.x
volumes = [analytic_volume, totcell, totvol, totint]
output("volumes")
for v in volumes:
    if np.abs(v - analytic_volume) > tolerance:
        print("volumes not correct")
        checksum += 1

totarea = 0.
for i in range(nodes.numNodes):
    totarea += np.sum(np.abs(scn[i]))
print("areas: ", totarea, analytic_surface)
if np.abs(totarea - analytic_surface) > tolerance:
    print("areas not correct")
    checksum += 1

#-------------------------------------------------------------------------------
# Check integrals
#-------------------------------------------------------------------------------
nPerhTest = 4.01
if (nx == 20) and (dimension == 1) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0):
    indi = 10
    indj = 11
    if useOverlap:
        indij = flatConnectivity.localToFlatOverlap(indi, indj)
    else:
        indij = flatConnectivity.localToFlat(indi, indj)
    if useOverlap:
        vals = [["vlK",  vlK[indi], 1.0],
                ["vlG",  vlG[indi].x, 0.0],
                ["vbKK",  vbKK[indi][indij], 0.9585478852898509],
                ["vbGK",  vbGK[indi][indij].x, -1.4389923272268597],
                ["vbKG",  vbKG[indi][indij].x, 1.4389923272268597],
                ["vbGdG",  vbGdG[indi][indij], 5.1346676217110305],
                ["vbGpG",  vbGpG[indi][indij].xx, 5.1346676217110305]]
    else:
        vals = [["vlK",  vlK[indi], 1.0],
                ["vlG",  vlG[indi].x, 0.0],
                ["vbKK",  vbKK[indi][indij], 1.21520485672],
                ["vbGK",  vbGK[indi][indij].x, -7.48643093476],
                ["vbKG",  vbKG[indi][indij].x, 7.48643093476],
                ["vbGdG",  vbGdG[indi][indij], -5.83078993373],
                ["vbGpG",  vbGpG[indi][indij].xx, -5.83078993373]]
    print("x: ", position(0, indi), position(0, indj))
    print("H: ", H(0, indi), H(0, indj))
    print("delta: ", delta[0])
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1

    indi = 0
    if useOverlap:
        indj = 2
        indij = flatConnectivity.localToFlatOverlap(indi, indj)
    else:
        indj = 1
        indij = flatConnectivity.localToFlat(indi, indj)
    numSurfaces = flatConnectivity.numSurfaces(indi)
    print("x: ", position(0, indi), position(0, indj))
    print("H: ", H(0, indi), H(0, indj))
    print("delta: ", 2*delta[0])
    if useOverlap:
        vals = [["slKn1",  slKn[indi][0].x, -0.7981466844744088],
                ["slKn2",  slKn[indj][0].x, -0.32244298020359935],
                ["slKKn",  sbKKn[indi][0 + numSurfaces * indij].x, -0.2573567955815503],
                ["vlK1",  vlK[indi], 0.581078921339981],
                ["vlK2",  vlK[indj], 0.9648661145429461],
                ["vlG1",  vlG[indi].x, -0.7981466844744085],
                ["vlG2",  vlG[indj].x, -0.32244298020359935],
                ["vbKK",  vbKK[indi][indij], 0.5297239342952187],
                ["vbGK",  vbGK[indi][indij].x, -0.6960555516515605],
                ["vbKG",  vbKG[indi][indij].x, 0.4386987560700106],
                ["vbGdG",  vbGdG[indi][indij], 0.691599920549981],
                ["vbGpG",  vbGpG[indi][indij].xx, 0.691599920549981]]
    else:
        vals = [["slKn1",  slKn[indi][0].x, -1.49474091258],
                ["slKn2",  slKn[indj][0].x, -0.697023258026],
                ["slKKn",  sbKKn[indi][0 + numSurfaces * indij].x, -1.04186918079],
                ["vlK1",  vlK[indi], 0.658577434997],
                ["vlK2",  vlK[indj], 0.934274660301],
                ["vlG1",  vlG[indi].x, -1.49474091258],
                ["vlG2",  vlG[indj].x, -0.697023258026],
                ["vbKK",  vbKK[indi][indij], 0.962387521061],
                ["vbGK",  vbGK[indi][indij].x, -2.26223953892],
                ["vbKG",  vbKG[indi][indij].x, 1.22037035812],
                ["vbGdG",  vbGdG[indi][indij], 4.06585331025],
                ["vbGpG",  vbGpG[indi][indij].xx, 4.06585331025]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1

if (nx == 10) and (ny == 10) and (dimension == 2) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0):
    indi = 5
    indj = 14
    print("xi/j: ", position(0, indi), position(0, indj))
    print("H: ", H(0, indi), H(0, indj))
    if useOverlap:
        indij = flatConnectivity.localToFlatOverlap(indi, indj)
    else:
        indij = flatConnectivity.localToFlat(indi, indj)
    normali = Vector(0.0, -1.0)
    inds = flatConnectivity.surfaceIndex(indi, normali)
    output("inds")
    numSurfaces = flatConnectivity.numSurfaces(indi)
    if useOverlap:
        vals = [["slKn1x",  slKn[indi][0].x, 0.0],
                ["slKn1y",  slKn[indi][0].y, -0.669533189156],
                ["slKn2x",  slKn[indj][0].x, 0.0],
                ["slKn2y",  slKn[indj][0].y, -0.384126497652],
                ["slKKnx",  sbKKn[indi][0 + numSurfaces * indij].x, 0.0],
                ["slKKny",  sbKKn[indi][0 + numSurfaces * indij].y, -0.121348978374],
                ["vlK1",  vlK[indi], 0.640055056],
                ["vlK2",  vlK[indj], 0.899817182704],
                ["vlG1x",  vlG[indi].x, 0.000422802024285],
                ["vlG1y",  vlG[indi].y, -0.669533225255],
                ["vlG2x",  vlG[indj].x, 0.0],
                ["vlG2y",  vlG[indj].y, -0.384126497454],
                ["vbKK",  vbKK[indi][indij], 0.197565771183],
                ["vbGKx",  vbGK[indi][indij].x, 0.144680555004],
                ["vbGKy",  vbGK[indi][indij].y, -0.194445006415],
                ["vbKGx",  vbKG[indi][indij].x, -0.144680555416],
                ["vbKGy",  vbKG[indi][indij].y, 0.0730960304324],
                ["vbGdG",  vbGdG[indi][indij], 0.432183556241],
                ["vbGpGxx",  vbGpG[indi][indij].xx, 0.2572459162768114],
                ["vbGpGxy",  vbGpG[indi][indij].xy, 0.052853608053438667],
                ["vbGpGyx",  vbGpG[indi][indij].yx, 0.14212702761020798],
                ["vbGpGyy",  vbGpG[indi][indij].yy, 0.17493765529791133]]
    else:
        vals = [["slKn1x",  slKn[indi][0].x, 0.0],
                ["slKn1y",  slKn[indi][0].y, -1.0962491626538495],
                ["slKn2x",  slKn[indj][0].x, 0.0],
                ["slKn2y",  slKn[indj][0].y, -0.059398685631869556],
                ["slKKnx",  sbKKn[indi][0 + numSurfaces * indij].x, 0.0],
                ["slKKny",  sbKKn[indi][0 + numSurfaces * indij].y, -0.03904333336926666],
                ["vlK1",  vlK[indi], 0.7597164007873803],
                ["vlK2",  vlK[indj], 0.99663461217954],
                ["vlG1x",  vlG[indi].x, 0.0],
                ["vlG1y",  vlG[indi].y, -1.0962502554318752],
                ["vlG2x",  vlG[indj].x, 0.0],
                ["vlG2y",  vlG[indj].y, -0.05940042563582875],
                ["vbKK",  vbKK[indi][indij], 0.36676470245842147],
                ["vbGKx",  vbGK[indi][indij].x, 1.0689592106857209],
                ["vbGKy",  vbGK[indi][indij].y, -1.08118303773501],
                ["vbKGx",  vbKG[indi][indij].x, -1.0689597769560295],
                ["vbKGy",  vbKG[indi][indij].y, 1.04213968246473],
                ["vbGdG",  vbGdG[indi][indij], -0.7581354835736813],
                ["vbGpGxx",  vbGpG[indi][indij].xx, -0.32659928253690645],
                ["vbGpGxy",  vbGpG[indi][indij].xy, 2.9109947113118233],
                ["vbGpGyx",  vbGpG[indi][indij].yx, 3.03962231782225],
                ["vbGpGyy",  vbGpG[indi][indij].yy, -0.4315362010367751]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1

if (nx == 5) and (ny == 5) and (nz == 5) and (dimension == 3) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0) and not useOverlap:
    indi = 30
    indj = 31
    print("xi/j: ", position(0, indi), position(0, indj))
    print("H: ", H(0, indi), H(0, indj))
    indij = flatConnectivity.localToFlat(indi, indj)
    normali1 = Vector(-1.0, 0.0, 0.0)
    normali2 = Vector(0.0, -1.0, 0.0)
    normali3 = Vector(0.0, 0.0, -1.0)
    inds1 = flatConnectivity.surfaceIndex(indi, normali1)
    inds2 = flatConnectivity.surfaceIndex(indi, normali2)
    inds3 = flatConnectivity.surfaceIndex(indi, normali3)
    numSurfaces = flatConnectivity.numSurfaces(indi)
    vals =  [["slKn1x",  slKn[indi][inds1].x, -0.5103947446431758],
             ["slKn2y",  slKn[indi][inds2].y, -0.06916083981314598],
             ["slKn3z",  slKn[indi][inds3].z, -0.06916083981314597],
             ["slKKn1x",  sbKKn[indi][inds1 + numSurfaces * indij].x, -0.007005827977086793],
             ["slKKn2y",  sbKKn[indi][inds2 + numSurfaces * indij].y, -0.0007475435998703369],
             ["slKKn3z",  sbKKn[indi][inds3 + numSurfaces * indij].z, -0.0007475435998703376],
             ["vlK1",  vlK[indi], 0.7159473522496711],
             ["vlK2",  vlK[indj], 0.9797160377111401],
             ["vlG1x",  vlG[indi].x, -0.5103944681812282],
             ["vlG1y",  vlG[indi].y, -0.06916049229770833],
             ["vlG1z",  vlG[indi].z, -0.06916049229770928],
             ["vbKK",  vbKK[indi][indij], 0.07601956302191366],
             ["vbGKx",  vbGK[indi][indij].x, -0.09426859868414547],
             ["vbGKy",  vbGK[indi][indij].y, -0.0002876532413641895],
             ["vbGKz",  vbGK[indi][indij].y, -0.0002876532413641895],
             ["vbKGx",  vbKG[indi][indij].x, 0.08726278511986846],
             ["vbKGy",  vbKG[indi][indij].y, -0.0004598944498155447],
             ["vbKGy",  vbKG[indi][indij].z, -0.0004598944498155505],
             ["vbGdG",  vbGdG[indi][indij], 0.22624447003774054],
             ["vbGpGxx",  vbGpG[indi][indij].xx, -2.0906918461268557e-05],
             ["vbGpGxy",  vbGpG[indi][indij].xy, 0.0007170740815179568],
             ["vbGpGxz",  vbGpG[indi][indij].xz, 0.0007170740815178934],
             ["vbGpGyx",  vbGpG[indi][indij].yx, -0.00042323317663459705],
             ["vbGpGyy",  vbGpG[indi][indij].yy, 0.11313268847810136],
             ["vbGpGyz",  vbGpG[indi][indij].yz, 1.7421837009401335e-07],
             ["vbGpGzx",  vbGpG[indi][indij].zx, -0.0004232331766347246],
             ["vbGpGzy",  vbGpG[indi][indij].zy, 1.742183700963299e-07],
             ["vbGpGzy",  vbGpG[indi][indij].zz, 0.11313268847810101]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1
            
#-------------------------------------------------------------------------------
# Check whether surface and volume integrals agree
#-------------------------------------------------------------------------------
print("surface-volume equivalence")
av_neighbors = 0.
av_surfaces = 0.
num_overlaps = 0
lin_err = [0., 0.]
bil_err = [0., 0.]
lin_sum = [0., 0.]
bil_sum = [0., 0.]
sv_time = time.time()
if checkAllSurfaceVolume:
    num_tested = nodes.numNodes
    ivals = list(range(nodes.numNodes))
else:
    num_tested = 20
    ivals = np.random.choice(list(range(nodes.numNodes)), size=num_tested)
    # num_tested = 2
    # ivals = [0, 1]

for i in ivals:
    if useOverlap:
        numElements = flatConnectivity.numOverlapNeighbors(i)
    else:
        numElements = flatConnectivity.numNeighbors(i)
    numSurfaces = flatConnectivity.numSurfaces(i)
    av_neighbors += numElements
    av_surfaces += numSurfaces
    sbKKni = sbKKn[i]
    sbKKn2i = sbKKn2[i]
    # print i, numElements, numSurfaces, position(0, i), flags(0, i), flatConnectivity.normal([0, i])
    # print sbKKni
    for j in range(numElements):
        t1 = sbKKn2i[j]
        t2 = Vector()
        for s in range(numSurfaces):
            t2 += sbKKni[s + numSurfaces * j]
        bil_err[0] += (t2 - t1).magnitude()
        bil_err[1] += (t2 - t1).magnitude2()
        bil_sum[0] += t2.magnitude()
        bil_sum[1] += t2.magnitude2()
        if outputAllSurfaceVolume:
            print("\tsbKKn\t{}\t{}\t{}\t{}\t{}".format(i, j, t1, t2, t2 - t1))
        num_overlaps += 1
    t1 = slKn2[i]
    t2 = Vector()
    slKni = slKn[i]
    for s in range(numSurfaces):
        t2 += slKni[s]
    lin_err[0] += (t2 - t1).magnitude()
    lin_err[1] += (t2 - t1).magnitude2()
    lin_sum[0] += t2.magnitude()
    lin_sum[1] += t2.magnitude2()
    if outputAllSurfaceVolume:
        print("\tslKn\ti\t{}\t{}\t{}".format(i, t1, t2, t2 - t1))
sv_time = time.time() - sv_time
av_neighbors /= num_tested
av_surfaces /= num_tested
for i in range(2):
    bil_err[i] /= (bil_sum[i] + 1.e-16)
    lin_err[i] /= (lin_sum[i] + 1.e-16)
bil_err[1] = np.sqrt(bil_err[1])
lin_err[1] = np.sqrt(lin_err[1])
output("av_neighbors")
output("av_surfaces")
for err in bil_err:
    if err > tolerance:
        checksum += 1
        print("bilinear error too high")
output("bil_err")
for err in lin_err:
    if err > tolerance:
        checksum += 1
        print("linear error too high")
output("lin_err")
output("sv_time")
    
#-------------------------------------------------------------------------------
# Output how many things have failed
#-------------------------------------------------------------------------------
output("checksum")
if checksum > 0:
    raise ValueError("too many errors")
