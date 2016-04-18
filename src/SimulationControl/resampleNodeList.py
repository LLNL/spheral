import mpi
import VoronoiDistributeNodes

#-------------------------------------------------------------------------------
# Resample to a new set of nodes represented by a generator.
#-------------------------------------------------------------------------------
def resampleNodeList(nodes,
                     generator,
                     W,
                     mask = None,
                     etaExclude = None,
                     boundaryConditions = []):

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

    # Check if we're doing a Solid or FluidNodeList.
    if isinstance(nodes, SolidNodeList):
        solid = True
    elif isinstance(nodes, FluidNodeList):
        solid = False
    else:
        raise RuntimeError, "Unknown NodeList type."

    # Build a temporary FluidNodeList we'll use to sample to.
    if solid:
        newnodes = makeSolidNodeList(name = "zznewnodes", 
                                     eos = nodes.eos,
                                     strength = nodes.strength,
                                     hmin = 1e-10,
                                     hmax = 1e10,
                                     xmin = -10*hmax*Vector.one,
                                     xmax =  10*hmax*Vector.one)
        if mask:
            masknodes = makeSolidNodeList(name = "zznewnodes", 
                                          eos = nodes.eos,
                                          strength = nodes.strength,
                                          hmin = 1e-10,
                                          hmax = 1e10,
                                          xmin = -10*hmax*Vector.one,
                                          xmax =  10*hmax*Vector.one)
    else:
        newnodes = makeFluidNodeList(name = "zznewnodes", 
                                     eos = nodes.eos,
                                     hmin = 1e-10,
                                     hmax = 1e10,
                                     xmin = -10*hmax*Vector.one,
                                     xmax =  10*hmax*Vector.one)
        if mask:
            masknodes = makeFluidNodeList(name = "zznewnodes", 
                                          eos = nodes.eos,
                                          hmin = 1e-10,
                                          hmax = 1e10,
                                          xmin = -10*hmax*Vector.one,
                                          xmax =  10*hmax*Vector.one)
    distributor((newnodes, generator))

    # Build the connectivity so we can figure out which nodes touch.
    db = DataBase()
    db.appendNodeList(nodes)
    db.appendNodeList(newnodes)
    nodes.neighbor().updateNodes()
    newnodes.neighbor().updateNodes()
    db.updateConnectivityMap(False)
    cm = db.connectivityMap()

    # If we're masking some points, things get complicated.  The mask nodes are going to persist to the new
    # nodes, and so we need to not overlay them.  We also want to remove any new nodes that overlap with the
    # mask nodes, since the masked ones are going to be copied to the new nodes in the end.
    if mask:

        # Copy the field values from the original masked nodes to the temporary mask set.
        

    oldkillNodes = vector_of_int()
    if mask is None:
        for i in xrange(nodes.numInternalNodes):
            oldnodes2kill.append(i)
    else:
        if etaExclude is None:
            etaExclude = nodes.neighbor().kernelExtent/2
        assert etaExclude > 0.0

        posi = nodes.positions()
        posj = newnodes.positions()
        H = nodes.Hfield()
        newnodes2kill = vector_of_int()
        for i in xrange(newnodes.numInternalNodes):
            fullconnectivity = cm.connectivityForNode(0, i)
            for j in fullconnectivity[1]:
                eta = (H[i]*(posi[i] - posj[j])).magnitude()
                if eta < etaExclude:
                    newnodes2kill.append(j)

        print "Removing %i nodes from new list due to overlap with masked nodes." % mpi.allreduce(len(newnodes2kill), mpi.SUM)
        newnodes.deleteNodes(newnodes2kill)

        # Update connectivity and such.
        newnodes.neighbor().updateNodes()
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()

    # Convert fields we're going to map to conserved values.  This is necessary 'cause the splat operation we're going
    # to use guarantees summing over the input and output field values gives the same value.
    mass = nodes.mass()
    momentum = VectorField(nodes.velocity())
    thermalenergy = ScalarField(nodes.specificThermalEnergy)
    for i in xrange(nodes.numNodes):
        momentum[i] *= mass[i]
        thermalenergy[i] *= mass[i]
    if isinstance(nodes, SolidNodeList):
        print "Copying solid fields."
        S = SymTensorField(nodes.deviatoricStress())
        ps = ScalarField(nodes.plasticStrain())
        D = SymTensorField(nodes.damage())
        for i in xrange(nodes.numNodes):
            S[i] *= mass[i]
            ps[i] *= mass[i]
            D[i] *= mass[i]

    # Map stuff from the old to new nodes.
    fls = FieldListSet()
    mass_fl = ScalarFieldList(Reference)
    momentum_fl = VectorFieldList(Reference)
    thermalenergy_fl = ScalarFieldList(Reference)
    fls.ScalarFieldLists.append(mass_fl)
    fls.VectorFieldLists.append(momentum_fl)
    fls.ScalarFieldLists.append(thermalenergy_fl)
    if isinstance(nodes, SolidNodeList):
        S_fl = SymTensorFieldList(Reference)
        ps_fl = ScalarFieldList(Reference)
        D_fl = SymTensorFieldList(Reference)
        fls.SymTensorFieldLists.append(S_fl)
        fls.ScalarFieldLists.append(ps_fl)
        fls.SymTensorFieldLists.append(D_fl)

    pos0_fl = VectorFieldList(Copy)
    mass0_fl = ScalarFieldList(Copy)
    H0_fl = SymTensorFieldList(Copy)
    pos0_fl.appendField(nodes.positions())
    mass0_fl.appendField(nodes.mass())
    H0_fl.appendField(nodes.Hfield())
    pos1_fl = VectorFieldList(Copy)
    mass1_fl = ScalarFieldList(Copy)
    H1_fl = SymTensorFieldList(Copy)
    pos1_fl.appendField(newnodes.positions())
    mass1_fl.appendField(newnodes.mass())
    H1_fl.appendField(newnodes.Hfield())
    bcs = vector_of_Boundary()
    for bc in boundaryConditions:
        bcs.append(bc)
    newfls = splatMultipleFieldsMash(fls,
                                     pos0_fl, mass0_fl, H0_fl, W,
                                     pos1_fl, mass1_fl, H1_fl, bcs)

    
    # Delete all the nodes we're no longer using from the old NodeList.
    killNodes = vector_of_int()
