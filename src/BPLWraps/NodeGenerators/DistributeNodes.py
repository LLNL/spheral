from math import *
import Spheral
import distributeNodesGeneric
import loadmpi
mpi, rank, numDomains = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Domain decompose using the sort and divide scheme (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(*listOfNodeTuples):
    if numDomains > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase1d,
                                                      Spheral.globalNodeIDs,
                                                      Spheral.SortAndDivideRedistributeNodes1d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase1d,
                                                      Spheral.globalNodeIDs1d,
                                                      None)

#-------------------------------------------------------------------------------
# Domain decompose using the sort and divide scheme (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    if numDomains > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase2d,
                                                      Spheral.globalNodeIDs,
                                                      Spheral.SortAndDivideRedistributeNodes2d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase2d,
                                                      Spheral.globalNodeIDs,
                                                      None)

#-------------------------------------------------------------------------------
# Domain decompose using the sort and divide scheme (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    if numDomains > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase3d,
                                                      Spheral.globalNodeIDs,
                                                      Spheral.SortAndDivideRedistributeNodes3d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase3d,
                                                      Spheral.globalNodeIDs,
                                                      None)

#-------------------------------------------------------------------------------
# Provide no-op versions of the distributer.
#-------------------------------------------------------------------------------
def nullDistributeNodes1d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase1d,
                                                  Spheral.globalNodeIDs,
                                                  None)

def nullDistributeNodes2d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase2d,
                                                  Spheral.globalNodeIDs,
                                                  None)

def nullDistributeNodes3d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase3d,
                                                  Spheral.globalNodeIDs,
                                                  None)

#-------------------------------------------------------------------------------
# Old method for handling 1-D distributions.  Deprecated now.
# Distribute sets of nodes evenly in the given ranges (1d)
# Input:  a list of tuples with (nodeList, (global) numNodes, (global) xRange)
#-------------------------------------------------------------------------------
def distributeNodesInRange1d(listOfNodeTuples,
                             nPerh = 2.01,
                             reverse = False):


    # Count how many nodes total we need to distribute.
    globalNumNodes = 0
    for (nodeList, n, rho, (x1, x2)) in listOfNodeTuples:
        assert n >= 0
        globalNumNodes += n

    # Determine the global node ID range this process will cover.
    nNodesPerDomain = globalNumNodes/numDomains
    minNodeID = rank*nNodesPerDomain
    if rank == numDomains - 1:
        maxNodeID = globalNumNodes
    else:
        maxNodeID = minNodeID + nNodesPerDomain
    
    # Set the size of each node list, and assign node positions.
    cumulativeN = 0
    for (nodeList, n, rho, (x0, x1)) in listOfNodeTuples:
        mi = (x1 - x0)*rho/max(1, n)
        minGlobalNodeListID = max(0, min(n, minNodeID - cumulativeN))
        maxGlobalNodeListID = max(0, min(n, maxNodeID - cumulativeN))
        nNodesThisDomain = maxGlobalNodeListID - minGlobalNodeListID
        assert(minGlobalNodeListID >= minGlobalNodeListID)
        nodeList.numInternalNodes = nNodesThisDomain
        cumulativeN += n
        dx = (x1 - x0)/n
        Hi = Spheral.SymTensor1d(1.0/(nPerh*dx))
        for globalNodeID in xrange(minGlobalNodeListID, maxGlobalNodeListID):
            localNodeID = globalNodeID - minGlobalNodeListID
            if reverse:
                nodeList.positions()[localNodeID].x = x1 - (globalNodeID + 0.5)*dx
            else:
                nodeList.positions()[localNodeID].x = x0 + (globalNodeID + 0.5)*dx
            nodeList.mass()[localNodeID] = mi
            nodeList.massDensity()[localNodeID] = rho
            nodeList.Hfield()[localNodeID] = Hi

        # Make sure the neighbor is up to date!
        nodeList.neighbor().updateNodes()
        assert nodeList.neighbor().valid()
        #nodeList.updateWeight()

    return

#-------------------------------------------------------------------------------
# Handle psuedo-1D spherical node distribution similarly to above.
# Input:  a list of tuples with (nodeList, (global) numNodes, (global) xRange)
#-------------------------------------------------------------------------------
def distributeNodesInSphericalRange3d(listOfNodeTuples,
                                      nPerh = 2.01,
                                      reverse = False):
    
    mpi, rank, numDomains = loadmpi.loadmpi()

    # Count how many nodes total we need to distribute.
    globalNumNodes = 0
    for (nodeList, n, rho, range) in listOfNodeTuples:
        assert n >= 0
        globalNumNodes += n

    # Determine the global node ID range this process will cover.
    nNodesPerDomain = globalNumNodes/numDomains
    minNodeID = rank*nNodesPerDomain
    if rank == numDomains - 1:
        maxNodeID = globalNumNodes
    else:
        maxNodeID = minNodeID + nNodesPerDomain
    
    # Set the size of each node list, and assign node positions, masses, and
    # H's.
    cumulativeN = 0
    for (nodeList, n, rho, (r0, r1)) in listOfNodeTuples:
        dr = (r1 - r0)/max(1, n)
        minGlobalNodeListID = max(0, min(n, minNodeID - cumulativeN))
        maxGlobalNodeListID = max(0, min(n, maxNodeID - cumulativeN))
        nNodesThisDomain = maxGlobalNodeListID - minGlobalNodeListID
        assert minGlobalNodeListID >= minGlobalNodeListID
        nodeList.numInternalNodes = nNodesThisDomain
        cumulativeN += n
        Hxx = 1.0/(nPerh*dr)
        Hi = Spheral.SymTensor3d(Hxx, 0.0, 0.0,
                                 0.0, Hxx, 0.0,
                                 0.0, 0.0, Hxx)

        vol = 4.0/3.0*pi*(r1**3 - r0**3)
        neff = vol/dr**3
        mi = vol*rho/neff

        for globalNodeID in xrange(minGlobalNodeListID, maxGlobalNodeListID):
            localNodeID = globalNodeID - minGlobalNodeListID
            if reverse:
                nodeList.positions()[localNodeID].x = r1 - (globalNodeID + 0.5)*dr
            else:
                nodeList.positions()[localNodeID].x = r0 + (globalNodeID + 0.5)*dr
            ri = nodeList.positions()[localNodeID].x

            nodeList.mass()[localNodeID] = mi
            nodeList.massDensity()[localNodeID] = rho
            nodeList.Hfield()[localNodeID] = Hi

        # Make sure the neighbor is up to date!
        nodeList.neighbor().updateNodes()
        assert nodeList.neighbor().valid()
        #nodeList.updateWeight()

    return

