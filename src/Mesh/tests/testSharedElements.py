import mpi

#----------------------------------------------------------------------
# Test that the shared nodes are consisten between domains.
#----------------------------------------------------------------------
def testSharedNodes(mesh):
    assert len(mesh.neighborDomains) == len(mesh.sharedNodes)

    # First check that everyone agrees about who is talking to who.
    for sendProc in xrange(mpi.procs):
        otherProcs = mpi.bcast(list(mesh.neighborDomains), root=sendProc)
        if mpi.rank != sendProc:
            assert (mpi.rank in otherProcs) == (sendProc in mesh.neighborDomains)

    # Build our set of global shared node IDs.
    globalIDs = mesh.globalMeshNodeIDs()
    globalSharedNodes = [[globalIDs[i] for i in localNodes] for localNodes in mesh.sharedNodes]
    assert len(globalSharedNodes) == len(mesh.neighborDomains)

    # Check that the shared nodes are consistent.
    sendRequests = []
    for (otherProc, ids) in zip(mesh.neighborDomains, globalSharedNodes):
        sendRequests.append(mpi.isend(ids, dest=otherProc))
    for (otherProc, ids) in zip(mesh.neighborDomains, globalSharedNodes):
        otherIDs = mpi.recv(source=otherProc)[0]
        assert ids == otherIDs

    return true
