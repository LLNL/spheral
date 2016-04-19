import mpi
import VoronoiDistributeNodes
import SolidSpheral

#...........................................................................
# A local helper method for copying data from one NodeList to another.
#...........................................................................
def copyNodeListFields(nodes0, nodes1, mask, solid):

    m0 = nodes0.mass()
    p0 = nodes0.positions()
    v0 = nodes0.velocity()
    H0 = nodes0.Hfield()
    r0 = nodes0.massDensity()
    e0 = nodes0.specificThermalEnergy()

    m1 = nodes1.mass()
    p1 = nodes1.positions()
    v1 = nodes1.velocity()
    H1 = nodes1.Hfield()
    r1 = nodes1.massDensity()
    e1 = nodes1.specificThermalEnergy()

    if solid:
        S0 = nodes0.deviatoricStress()
        ps0 = nodes0.plasticStrain()
        psr0 = nodes0.plasticStrainRate()
        D0 =  nodes0.damage()
        
        S1 = nodes1.deviatoricStress()
        ps1 = nodes1.plasticStrain()
        psr1 = nodes1.plasticStrainRate()
        D1 =  nodes1.damage()

    j = 0
    for i in xrange(nodes0.numInternalNodes):
        if mask[i] == 1:
            assert j < nodes1.numInternalNodes
            m1[j] = m0[i]
            p1[j] = p0[i]
            v1[j] = v0[i]
            H1[j] = H0[i]
            r1[j] = r0[i]
            e1[j] = e0[i]
            if solid:
                S1[j] = S0[i]
                ps1[j] = ps0[i]
                psr1[j] = psr0[i]
                D1[j] = D0[i]
            j += 1
    return

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
    if isinstance(nodes, SolidSpheral.NodeList1d):
        ndim = 1
    elif isinstance(nodes, SolidSpheral.NodeList2d):
        ndim = 2
    elif isinstance(nodes, SolidSpheral.NodeList3d):
        ndim = 3
    else:
        raise ValueError, "Unknown thing %s handed in: expected a NodeList" % nodes
    ndim0 = ndim
    exec "from SolidSpheral%id import *" % ndim   # Load the aliases for our dimensionality
    ndim = ndim0
    exec "from VoronoiDistributeNodes import distributeNodes%id as distributor" % ndim

    # Check if we're doing a Solid or FluidNodeList.
    if isinstance(nodes, SolidNodeList):
        solid = True
        NLF = makeSolidNodeList
    elif isinstance(nodes, FluidNodeList):
        solid = False
        NLF = makeFluidNodeList
    else:
        raise RuntimeError, "Unknown NodeList type."

    # Build a temporary NodeList we'll use to sample to.
    newnodes = NLF(name = "zza_newnodes", 
                   eos = nodes.eos,
                   hmin = 1e-10,
                   hmax = 1e10,
                   xmin = -10*nodes.hmax*Vector.one,
                   xmax =  10*nodes.hmax*Vector.one)
    if mask:
        masknodes = NLF(name = "zzz_masknodes", 
                        eos = nodes.eos,
                        hmin = 1e-10,
                        hmax = 1e10,
                        xmin = -10*nodes.hmax*Vector.one,
                        xmax =  10*nodes.hmax*Vector.one)
    distributor((newnodes, generator))

    # If we're masking some points, things get complicated.  The mask nodes are going to persist to the new
    # nodes, and so we need to not overlay them.  We also want to remove any new nodes that overlap with the
    # mask nodes, since the masked ones are going to be copied to the new nodes in the end.
    nmask = 0
    if mask:

        # Copy the field values from the original masked nodes to the temporary mask set.
        nmask = mask.localSum()
        print "Copying %i masked nodes from the original NodeList." % nmask
        masknodes.numInternalNodes = nmask
        copyNodeListFields(nodes, masknodes, mask, solid)

        # Remove the mask nodes from the starting NodeList.
        nodes2kill = vector_of_int()
        for i in xrange(nodes.numInternalNodes):
            if mask[i] == 1:
                nodes2kill.append(i)
        assert nodes2kill.size() == nmask
        nodes.deleteNodes(nodes2kill)

        # Now we need to remove any nodes from the target set that overlap with the mask nodes.
        db = DataBase()
        db.appendNodeList(newnodes)
        db.appendNodeList(masknodes)
        newnodes.neighbor().updateNodes()
        masknodes.neighbor().updateNodes()
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()

        if etaExclude is None:
            etaExclude = nodes.neighbor().kernelExtent/2.0
        assert etaExclude > 0.0

        posmask = masknodes.positions()
        Hmask = masknodes.Hfield()
        posnew = newnodes.positions()
        Hnew = newNodes.Hfield()
        nodes2kill = vector_of_int()
        for i in xrange(masknodes.numInternalNodes):
            fullconnectivity = cm.connectivityForNode(1, i)
            for j in fullconnectivity[0]:
                eta = min((Hmask[i]*(posmask[i] - posnew[j])).magnitude(),
                          (Hnew[j]*(posmask[i] - posnew[j])).magnitude())
                if eta < etaExclude:
                    nodes2kill.append(j)

        print "Removing %i nodes from new list due to overlap with masked nodes." % mpi.allreduce(len(nodes2kill), mpi.SUM)
        newnodes.deleteNodes(nodes2kill)

    # Build the connectivity so we can do the overlay.
    db = DataBase()
    db.appendNodeList(nodes)
    db.appendNodeList(newnodes)
    nodes.neighbor().updateNodes()
    newnodes.neighbor().updateNodes()
    # db.updateConnectivityMap(False)
    # cm = db.connectivityMap()

    # Convert fields we're going to map to conserved values.  This is necessary 'cause the splat operation we're going
    # to use guarantees summing over the input and output field values gives the same value.
    mass = nodes.mass()
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    momentum = VectorField(vel)
    thermalenergy = ScalarField(eps)
    for i in xrange(nodes.numNodes):
        momentum[i] *= mass[i]
        thermalenergy[i] *= mass[i]
    if solid:
        S = nodes.deviatoricStress()
        ps = nodes.plasticStrain()
        D = nodes.damage()
        mS = SymTensorField(S)
        mps = ScalarField(ps)
        mD = SymTensorField(D)
        for i in xrange(nodes.numNodes):
            mS[i] *= mass[i]
            mps[i] *= mass[i]
            mD[i] *= mass[i]

    # Map stuff from the old to new nodes.
    fls = FieldListSet()
    mass_fl = ScalarFieldList()
    momentum_fl = VectorFieldList()
    thermalenergy_fl = ScalarFieldList()
    mass_fl.appendField(mass)
    momentum_fl.appendField(momentum)
    thermalenergy_fl.appendField(thermalenergy)
    mass_fl.copyFields()
    momentum_fl.copyFields()
    thermalenergy_fl.copyFields()
    fls.ScalarFieldLists.append(mass_fl)
    fls.VectorFieldLists.append(momentum_fl)
    fls.ScalarFieldLists.append(thermalenergy_fl)
    if solid:
        S_fl = SymTensorFieldList()
        ps_fl = ScalarFieldList()
        D_fl = SymTensorFieldList()
        S_fl.appendField(mS)
        ps_fl.appendField(mps)
        D_fl.appendField(mD)
        S_fl.copyFields()
        ps_fl.copyFields()
        D_fl.copyFields()
        fls.SymTensorFieldLists.append(S_fl)
        fls.ScalarFieldLists.append(ps_fl)
        fls.SymTensorFieldLists.append(D_fl)

    pos0_fl = VectorFieldList()
    mass0_fl = ScalarFieldList()
    H0_fl = SymTensorFieldList()
    pos0_fl.appendField(nodes.positions())
    mass0_fl.appendField(nodes.mass())
    H0_fl.appendField(nodes.Hfield())
    pos1_fl = VectorFieldList()
    mass1_fl = ScalarFieldList()
    H1_fl = SymTensorFieldList()
    pos1_fl.appendField(newnodes.positions())
    mass1_fl.appendField(newnodes.mass())
    H1_fl.appendField(newnodes.Hfield())
    pos0_fl.copyFields()
    mass0_fl.copyFields()
    H0_fl.copyFields()
    pos1_fl.copyFields()
    mass1_fl.copyFields()
    H1_fl.copyFields()
    bcs = vector_of_Boundary()
    for bc in boundaryConditions:
        bcs.append(bc)
    print "Splatting fields..."
    # newfls = sampleMultipleScalarFieldsMash(fls,
    #                                  pos0_fl, mass0_fl, H0_fl, W,
    #                                  pos1_fl, mass1_fl, H1_fl)
    newfls = splatMultipleFieldsMash(fls,
                                     pos0_fl, mass0_fl, H0_fl, W,
                                     pos1_fl, mass1_fl, H1_fl)
    print "Done splatting."

    # Denormalize the mapped values and fill them in as new values for the nodes.
    nodes.numInternalNodes = nmask + newnodes.numInternalNodes
    print nodes.numInternalNodes, newnodes.numInternalNodes
    mass1 = newfls.ScalarFieldLists[0][0]
    momentum1 = newfls.VectorFieldLists[0][0]
    thermalenergy1 = newfls.ScalarFieldLists[1][0]
    for i in xrange(newnodes.numInternalNodes):
        j = nmask + i
        assert mass1[i] > 0.0
        print j, mass[j], mass1[i]
        mass[j] = mass1[i]
        vel[j] = momentum1[i]/mass1[i]
        eps[j] = thermalenergy1[i]/mass1[i]
    if solid:
        mS1 = newfls.SymTensorFieldLists[0][0]
        mps1 = newfls.ScalarFieldLists[2][0]
        mD1 = newfls.SymTensorFieldLists[1][0]
        for i in xrange(newnodes.numInternalNodes):
            j = nmask + i
            assert mass1[i] > 0.0
            S[j] = mS1[i]/mass1[i]
            ps[j] = msp1[i]/mass1[i]
            D[j] = mD1[i]/mass1[i]

    # Insert any masked nodes, and we're done.
    if mask:
        newmask = [1]*nmask + [0]*nodes.numInternalNodes
        copyNodeListFields(masknodes, nodes, newmask, solid)

    # Whew!
    print "Finished resampling nodes: final node count %i." % mpi.allreduce(nodes.numInternalNodes, mpi.SUM)
    return
