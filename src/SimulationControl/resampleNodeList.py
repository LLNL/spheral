import mpi
import VoronoiDistributeNodes

#-------------------------------------------------------------------------------
# Resample to a new set of nodes represented by a generator.
#-------------------------------------------------------------------------------
def resampleNodeList(nodes,
                     generator,
                     mask = None,
                     etaExclude = None):

    # Check our dimensionality
    if isinstance(nodes, NodeList1d):
        ndim = 1
    elif isinstance(nodes, NodeList2d):
        ndim = 2
    elif isinstance(nodes, NodeList3d):
        ndim = 3
    else:
        raise ValueError, "Unknown thing %s handed in: expected a NodeList" % nodes
    exec("from SolidSpheral%id import *" % ndim)   # Load the aliases for our dimensionality
    exec("from VoronoiDistributeNodes import distributeNodes%id as distributor" % ndim)

    # Build a temporary FluidNodeList we'll use to sample to.
    tmpnodes = makeSolidNodeList(name = "zztmpnodes", 
                                 eos = nodes.eos,
                                 hmin = 1e-10,
                                 hmax = 1e10,
                                 xmin = -10*hmax*Vector.one,
                                 xmax =  10*hmax*Vector.one)
    distributor((tmpnodes, generator))

    # Build the connectivity so we can figure out which nodes touch.
    db = DataBase()
    db.appendNodeList(nodes)
    db.appendNodeList(tmpnodes)
    nodes.neighbor().updateNodes()
    tmpnodes.neighbor().updateNodes()
    db.updateConnectivityMap(False)
    cm = db.connectivityMap()

    # Look for any new nodes we need to kill based on the mask.
    if not mask is None:
        if etaExclude is None:
            etaExclude = nodes.neighbor().kernelExtent/2
        assert etaExclude > 0.0

        posi = nodes.positions()
        posj = tmpnodes.positions()
        H = nodes.Hfield()
        nodes2kill = vector_of_int()
        for i in xrange(tmpnodes.numInternalNodes):
            fullconnectivity = cm.connectivityForNode(0, i)
            for j in fullconnectivity[1]:
                eta = (H[i]*(posi[i] - posj[j])).magnitude()
                if eta < etaExclude:
                    nodes2kill.append(j)

        print "Removing %i nodes from new list due to overlap with masked nodes." % mpi.allreduce(len(nodes2kill), mpi.SUM)
        tmpnodes.deleteNodes(nodes2kill)

        # Update connectivity and such.
        tmpnodes.neighbor().updateNodes()
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()

    # Convert fields we're going to map to conserved values.
    

    # Map stuff from the old to new nodes.
