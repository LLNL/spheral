import mpi
from Spheral2d import *

#----------------------------------------------------------------------
# Test that the shared nodes are consisten between domains.
#----------------------------------------------------------------------
def testSharedNodes(mesh):
    assert len(mesh.neighborDomains) == len(mesh.sharedNodes)

    # First check that everyone agrees about who is talking to who.
    myNeighborDomains = list(mesh.neighborDomains)
    for sendProc in xrange(mpi.procs):
        otherProcs = mpi.bcast(myNeighborDomains, root=sendProc)
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

    # Check that all shared nodes have been found.
    localSharedNodes = [[i for i in localNodes] for localNodes in mesh.sharedNodes]
    positions = vector_of_Vector()
    for i in xrange(mesh.numNodes):
        positions.append(mesh.node(i).position())
    xmin, xmax = Vector(), Vector()
    boundingBox(positions, xmin, xmax)
    xmin = Vector(mpi.allreduce(xmin.x, mpi.MIN), mpi.allreduce(xmin.y, mpi.MIN))
    xmax = Vector(mpi.allreduce(xmax.x, mpi.MAX), mpi.allreduce(xmax.y, mpi.MAX))
    boxInv = Vector(1.0/(xmax.x - xmin.x),
                    1.0/(xmax.y - xmin.y))
    nodeHashes = [hashPosition(mesh.node(i).position(), xmin, xmax, boxInv) for i in xrange(mesh.numNodes)]
    nodeHashes2ID = {}
    for i in xrange(len(nodeHashes)):
        nodeHashes2ID[nodeHashes[i]] = i
    for sendProc in xrange(mpi.procs):
        otherNodeHashes = mpi.bcast(nodeHashes, root=sendProc)
        if sendProc != mpi.rank:
            for hashi in otherNodeHashes:
                if hashi in nodeHashes:
                    assert sendProc in myNeighborDomains
                    idomain = myNeighborDomains.index(sendProc)
                    i = nodeHashes2ID[hashi]
                    assert i in localSharedNodes[idomain]

    # Same for faces.
    localSharedFaces = [[i for i in localFaces] for localFaces in mesh.sharedFaces]
    positions = vector_of_Vector()
    for i in xrange(mesh.numFaces):
        positions.append(mesh.face(i).position())
    faceHashes = [hashPosition(mesh.face(i).position(), xmin, xmax, boxInv) for i in xrange(mesh.numFaces)]
    faceHashes2ID = {}
    for i in xrange(len(faceHashes)):
        faceHashes2ID[faceHashes[i]] = i
    for sendProc in xrange(mpi.procs):
        otherFaceHashes = mpi.bcast(faceHashes, root=sendProc)
        if sendProc != mpi.rank:
            for hashi in otherFaceHashes:
                if hashi in faceHashes:
                    assert sendProc in myNeighborDomains
                    idomain = myNeighborDomains.index(sendProc)
                    i = faceHashes2ID[hashi]
                    assert i in localSharedFaces[idomain]

    return True
