import Spheral
import mpi

#-------------------------------------------------------------------------------
# Domain decompose using some specified domain partitioner (generic method).
#-------------------------------------------------------------------------------
def distributeNodesGeneric(listOfNodeTuples,
                           DataBaseType,
                           globalNodeIDs,
                           RedistributeNodesType):

    # We'll build the NodeLists into a DataBase.
    db = DataBaseType()

    # Assign nodes to domains by globalNodeID as a first cut.
    kernelExtent = 0.0
    numNodesPerProcess = [0]
    totalNumGlobalNodes = 0
    for (nodes, generator) in listOfNodeTuples:
        nglobal = generator.globalNumNodes()
        nlocal = generator.localNumNodes()
        print "distributeNodesGeneric: working on %s, (local, global) number nodes %i %i" % (nodes.name, nlocal, nglobal)
        numNodesPerProcess[0] += nlocal
        totalNumGlobalNodes += nglobal

        # Find the maximum kernel extent for all NodeLists.
        kernelExtent = max(kernelExtent, nodes.neighbor().kernelExtent)
        one = eval("Spheral.SymTensor%id.one" % db.nDim)
        hminInv = 1.0/nodes.hmin
        hmaxInv = 1.0/nodes.hmax

        # We start with the initial crappy distribution used in the generator.
        nodes.numGhostNodes = 0
        nodes.numInternalNodes = nlocal
        assert mpi.allreduce(nodes.numInternalNodes, mpi.SUM) == nglobal
        print "  distributeNodesGeneric: performing initial crappy distribution."
        r = nodes.positions()
        m = nodes.mass()
        vel = nodes.velocity()
        rho = nodes.massDensity()
        H = nodes.Hfield()
        for i in xrange(nlocal):
            r[i] = generator.localPosition(i)
            m[i] = generator.localMass(i)
            vel[i] = generator.localVelocity(i)
            rho[i] = generator.localMassDensity(i)
            H[i] = generator.localHtensor(i)
        H.applyScalarMin(hmaxInv)
        H.applyScalarMax(hminInv)

        # Put this NodeList into the DataBase.
        db.appendNodeList(nodes)
        print "  distributeNodesGeneric: %s initially finished" % nodes.name

    # # Update Neighbor information.
    # exec("Spheral.Neighbor%id.setBoundingBox()" % db.nDim)
    # for (nodes, generator) in listOfNodeTuples:
    #     nodes.neighbor().updateNodes()
    #     if (isinstance(nodes, Spheral.FluidNodeList1d) or
    #         isinstance(nodes, Spheral.FluidNodeList2d) or
    #         isinstance(nodes, Spheral.FluidNodeList3d)):
    #         nodes.updateWeight()

    # Report the initial breakdown.
    numNodesPerProcess = mpi.allreduce(numNodesPerProcess, mpi.SUM)
    print "(min, max, avg) nodes per process initially:  ", min(numNodesPerProcess), max(numNodesPerProcess), sum(numNodesPerProcess)/len(numNodesPerProcess)
    print "Total number of nodes: ", totalNumGlobalNodes

    # Now have the Redistributer repartition the nodes into something sensible.  Note this
    # automatically redistributes the globalNodeListID fields as well.
    print "distributeNodesGeneric: calling for redistribution."
    if RedistributeNodesType:
        repartition = RedistributeNodesType(kernelExtent)
        repartition.redistributeNodes(db)
    print "distributeNodesGeneric: redistribution done."

    # Update the neighboring info.
    #exec("Spheral.Neighbor%id.setBoundingBox()" % db.nDim)
    for nodes in db.nodeLists():
        nodes.neighbor().updateNodes()

    # Make sure we finished with the correct numbers of nodes!
    totalCheck = mpi.allreduce(sum([nodes.numInternalNodes for nodes in db.nodeLists()]), mpi.SUM)
    assert totalCheck == totalNumGlobalNodes
