from LLNLSpheral import *
from SpheralTestUtilities import *
import os, shutil, time

title("distributed connectivity")

commandLine(
    # Discretization options
    dimension = 1,
    nx = 16,
    x0 = 0.0,
    x1 = 1.0,
    nPerh = 4.01,

    # Optionally randomize the node positions
    randomizeNodes = False,
    ranfrac = 0.2,

    # Testing options
    testOverlap = False,
    printErrors = False,
    standardGlobalIDs = False,

    # Base directory for storing connectivity
    baseTestDir = "data-distributed")

if dimension == 1:
    from LLNLSpheral1d import *
elif dimension == 2:
    from LLNLSpheral2d import *
else:
    from LLNLSpheral3d import *

#-------------------------------------------------------------------------------
# Set up the output directories
#-------------------------------------------------------------------------------
dataDir = os.path.join(baseTestDir,
                       "dim={}".format(dimension),
                       "nx={}".format(nx),
                       "x0,x1={},{}".format(x0, x1),
                       "nPerh={}".format(nPerh),
                       "overlap={}".format(testOverlap),
                       "glIDs={}".format(standardGlobalIDs),
                       "rand={}".format(randomizeNodes if not randomizeNodes else "{},{}".format(randomizeNodes, ranfrac)))
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
        
#-------------------------------------------------------------------------------
# Interpolation kernels
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")
output("kernelExtent")

delta = (x1 - x0) / nx
hmid = delta * nPerh * kernelExtent
hmin = hmid * 1.e-3
hmax = hmid * 1.e3

#-------------------------------------------------------------------------------
# Material properties
#-------------------------------------------------------------------------------
units = MKS()
gamma = 5./3.
mu = 1.0
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Make the NodeList
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
# Seed the nodes
#-------------------------------------------------------------------------------
rho0 = 1.0
if dimension == 1:
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (x0, x1))],
                             nPerh = nPerh)
elif dimension == 2:
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(distributionType="lattice",
                                           nRadial = nx, nTheta = nx,
                                           xmin = (x0, x0),
                                           xmax = (x1, x1),
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
                                           n1 = nx, n2 = nx, n3 = nx,
                                           xmin = (x0, x0, x0),
                                           xmax = (x1, x1, x1),
                                           rho=rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes, generator))

#-------------------------------------------------------------------------------
# Randomize the node positions
#-------------------------------------------------------------------------------
if randomizeNodes:
    import random
    seed = 2
    rangen = random.Random()
    rangen.seed(seed)
    delta = (x1 - x0) / nx
    pos = nodes.positions()
    for i in xrange(nodes.numInternalNodes):
        if dimension == 1:
            pos[i].x += ranfrac * delta * rangen.uniform(-1.0, 1.0)
        elif dimension == 2:
            pos[i].x += ranfrac * delta * rangen.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * delta * rangen.uniform(-1.0, 1.0)
        elif dimension == 3:
            pos[i].x += ranfrac * delta * rangen.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * delta * rangen.uniform(-1.0, 1.0)
            pos[i].z += ranfrac * delta * rangen.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Make the DataBase
#-------------------------------------------------------------------------------
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")
output("dataBase.numNodeLists")
output("dataBase.numFluidNodeLists")
output("dataBase.numNodes")
output("dataBase.numInternalNodes")

#-------------------------------------------------------------------------------
# Get some fields from DataBase
#-------------------------------------------------------------------------------
mass = dataBase.globalMass
position = dataBase.globalPosition
h = dataBase.globalHfield
fieldLists = [mass, position, h]

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
dataBase.updateConnectivityMap(testOverlap, testOverlap)

#-------------------------------------------------------------------------------
# Create the distributed boundary and compute the connectivity
#-------------------------------------------------------------------------------
domainbc = TreeDistributedBoundary.instance()
domainbc.setAllGhostNodes(dataBase)
for fl in fieldLists:
    domainbc.applyFieldListGhostBoundary(fl)
domainbc.finalizeGhostBoundary()
for NN in dataBase.nodeLists():
    NN.neighbor().updateNodes()
dataBase.updateConnectivityMap(testOverlap, testOverlap)
connectivity = dataBase.connectivityMap(testOverlap, testOverlap)
output("connectivity")
output("dataBase.numNodes")
output("dataBase.numInternalNodes")

#-------------------------------------------------------------------------------
# Compute the global node IDs
#-------------------------------------------------------------------------------
if standardGlobalIDs:
    if dimension > 1:
        print "warning: test will likely fail, as global IDs depends on mpi.procs"
    globalIndicesFL = globalNodeIDsAll(dataBase)
else:
    globalIndicesFL = mortonOrderIndices(dataBase)
    # globalIndicesFL = peanoHilbertOrderIndices(dataBase) # not working
for f in globalIndicesFL:
    domainbc.applyGhostBoundary(f)
domainbc.finalizeGhostBoundary()
output("globalIndicesFL")

#-------------------------------------------------------------------------------
# Get connectivity in terms of global indices
#-------------------------------------------------------------------------------
globalIndices = []
globalConnectivity = []
for ni in range(dataBase.numNodeLists):
    for i in range(connectivity.numNodes(ni)):
        globali = globalIndicesFL(ni, i)
        globalConnectivityi = []
        connectivityi = connectivity.overlapConnectivityForNode(ni, i) if testOverlap else connectivity.connectivityForNode(ni, i)
        for nj in range(dataBase.numNodeLists):
            for j in connectivityi[nj]:
                globalj = globalIndicesFL(nj, j)
                globalConnectivityi.append(globalj)
        globalIndices.append(globali)
        globalConnectivity.append(set(globalConnectivityi))
if mpi.procs > 1:
    globalIndices = mpi.allreduce(globalIndices)
    globalConnectivity = mpi.allreduce(globalConnectivity)
globalData = dict(zip(globalIndices, globalConnectivity))

#-------------------------------------------------------------------------------
# Save connectivity if serial, load regardless
#-------------------------------------------------------------------------------
import pickle
dataFilename = os.path.join(dataDir, "connectivity")
if mpi.procs == 1:
    with open(dataFilename, 'w') as f:
        pickle.dump(globalData, f)

if os.path.exists(dataFilename):
    with open(dataFilename, 'r') as f:
        serialData = pickle.load(f)
else:
    raise IOError, "need to run in serial first"

#-------------------------------------------------------------------------------
# Check connectivity
#-------------------------------------------------------------------------------e
if mpi.rank == 0:
    # Go through serial connectivity and make sure parallel connectivity is the same
    numFailures = 0
    for i, cs in serialData.items():
        if i in globalData:
            cp = globalData[i]
            # The parallel connectivity contains more points than the serial
            # (is this expected?), so we check whether the parallel connectivity
            # contains the complete serial connectivity
            if not cp.issuperset(cs):
                if printErrors:
                    diff = cp.difference(cs)
                    print "diff for point {}: {}".format(i, diff)
                numFailures += 1
        else:
            numFailures += 1
            print "point {} not found".format(i)
    if numFailures > 0:
        raise ValueError, "num points with connectivity differences: {}".format(numFailures)
    else:
        print "test passed"
        
    

    
