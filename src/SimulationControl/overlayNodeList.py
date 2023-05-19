import mpi
import time
import VoronoiDistributeNodes
import SolidSpheral

#...........................................................................
# A local helper method for copying data from one NodeList to another.
#...........................................................................
def copyNodeListFields(nodes0, nodes1, solid):

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

    return

#-------------------------------------------------------------------------------
# Resample to a new set of nodes represented by a generator.
#-------------------------------------------------------------------------------
def overlayNodeList(nodes,
                    generator,
                    removeUnusedNodes = True):

    # Check our dimensionality
    if isinstance(nodes, SolidSpheral.NodeList1d):
        ndim = 1
    elif isinstance(nodes, SolidSpheral.NodeList2d):
        ndim = 2
    elif isinstance(nodes, SolidSpheral.NodeList3d):
        ndim = 3
    else:
        raise ValueError("Unknown thing %s handed in: expected a NodeList" % nodes)
    ndim0 = ndim
    exec("from SolidSpheral%id import *" % ndim)   # Load the aliases for our dimensionality
    ndim = ndim0
    exec("from VoronoiDistributeNodes import distributeNodes%id as distributor" % ndim)

    # Clear out any initial ghost nodes.
    nodes.numGhostNodes = 0

    # Check if we're doing a Solid or FluidNodeList.
    if isinstance(nodes, SolidNodeList):
        solid = True
        NLF = makeSolidNodeList
    elif isinstance(nodes, FluidNodeList):
        solid = False
        NLF = makeFluidNodeList
    else:
        raise RuntimeError("Unknown NodeList type.")

    # Check how to set the new neighbor info.
    boundaries = vector_of_Boundary()
    if isinstance(nodes._neighbor, NestedGridNeighbor):
        topGridSize = nodes._neighbor.topGridSize
        xmin = Vector.zero
        xmax = Vector.one * topGridSize
        NeighborType = NestedGridNeighbor
        if mpi.procs > 1:
            dbc = NestedGridDistributedBoundary.instance()
            boundaries.append(dbc)
    elif isinstance(nodes._neighbor, TreeNeighbor):
        xmin = nodes._neighbor.xmin
        xmax = nodes._neighbor.xmax
        topGridSize = (xmax - xmin).maxAbsElement()
        NeighborType = TreeNeighbor
        if mpi.procs > 1:
            dbc = BoundingVolumeDistributedBoundary.instance()
            boundaries.append(dbc)
            #raise RuntimeError, "Need a parallel policy for TreeNeighbor."
    else:
        raise RuntimeError("Unknown Neighbor type.")

    # Build a temporary NodeList we'll use to sample to.
    newnodes = NLF(name = "zza_newnodes", 
                   eos = nodes.eos,
                   hmin = 1e-10,
                   hmax = 1e10,
                   NeighborType = NeighborType,
                   topGridCellSize = topGridSize,
                   xmin = xmin,
                   xmax = xmax)
    distributor((newnodes, generator))

    # Convert fields we're going to map to conserved values.  This is necessary 'cause the splat operation we're going
    # to use guarantees summing over the input and output field values gives the same value.
    pos0 = nodes.positions()
    H0 = nodes.Hfield()
    mass0 = nodes.mass()
    rho0 = nodes.massDensity()
    vol0 = ScalarField(nodes.mass())
    vel0 = nodes.velocity()
    eps0 = nodes.specificThermalEnergy()
    momentum0 = VectorField(vel0)
    thermalenergy0 = ScalarField(eps0)
    for i in range(nodes.numNodes):
        vol0[i] /= rho0[i] + 1.0e-30
        momentum0[i] *= mass0[i]
        thermalenergy0[i] *= mass0[i]
    if solid:
        S0 = nodes.deviatoricStress()
        ps0 = nodes.plasticStrain()
        D0 = nodes.damage()
        mS0 = SymTensorField(S0)
        mps0 = ScalarField(ps0)
        mD0 = SymTensorField(D0)
        for i in range(nodes.numNodes):
            mS0[i] *= mass0[i]
            mps0[i] *= mass0[i]
            mD0[i] *= mass0[i]
    pos1 = newnodes.positions()
    H1 = newnodes.Hfield()
    mass1 = newnodes.mass()
    rho1 = newnodes.massDensity()
    vol1 = ScalarField(newnodes.mass())
    vel1 = newnodes.velocity()
    eps1 = newnodes.specificThermalEnergy()
    momentum1 = VectorField(vel1)
    thermalenergy1 = ScalarField(eps1)
    if solid:
        S1 = newnodes.deviatoricStress()
        ps1 = newnodes.plasticStrain()
        D1 = newnodes.damage()
        mS1 = SymTensorField(S1)
        mps1 = ScalarField(ps1)
        mD1 = SymTensorField(D1)

    # Map stuff from the old to new nodes.
    donorScalarFields = vector_of_ScalarFieldPtr()
    donorVectorFields = vector_of_VectorFieldPtr()
    donorTensorFields = vector_of_TensorFieldPtr()
    donorSymTensorFields = vector_of_SymTensorFieldPtr()
    acceptorScalarFields = vector_of_ScalarFieldPtr()
    acceptorVectorFields = vector_of_VectorFieldPtr()
    acceptorTensorFields = vector_of_TensorFieldPtr()
    acceptorSymTensorFields = vector_of_SymTensorFieldPtr()
    for f in (mass0, vol0, thermalenergy0):
        donorScalarFields.append(f)
    donorVectorFields.append(momentum0)
    for f in (mass1, vol1, thermalenergy1):
        acceptorScalarFields.append(f)
    acceptorVectorFields.append(momentum1)
    if solid:
        donorSymTensorFields.append(mS0)
        donorScalarFields.append(mps0)
        donorSymTensorFields.append(mD0)
        acceptorSymTensorFields.append(mS1)
        acceptorScalarFields.append(mps1)
        acceptorSymTensorFields.append(mD1)
    print("Splatting fields...")
    t0 = time.time()
    overlayRemapFields(boundaries,
                       donorScalarFields, donorVectorFields, donorTensorFields, donorSymTensorFields,
                       acceptorScalarFields, acceptorVectorFields, acceptorTensorFields, acceptorSymTensorFields)
    print("Done splatting, required %g seconds." % (time.time() - t0))

    # Denormalize the mapped values and fill them in as new values for the nodes.
    nodes.numInternalNodes = newnodes.numInternalNodes
    for i in range(newnodes.numInternalNodes):
        j = i
        pos0[j] = pos1[i]
        H0[j] = H1[i]
        if mass1[i] > 0.0:
            assert vol1[i] > 0.0
            mass0[j] = mass1[i]
            rho0[j] = mass1[i]/vol1[i]
            vel0[j] = momentum1[i]/mass1[i]
            eps0[j] = thermalenergy1[i]/mass1[i]
    if solid:
        for i in range(newnodes.numInternalNodes):
            j = i
            if mass1[i] > 0.0:
                S0[j] = mS1[i]/mass1[i]
                ps0[j] = mps1[i]/mass1[i]
                D0[j] = mD1[i]/mass1[i]

    # Look for any nodes that didn't get any information in the new set and delete them.
    if removeUnusedNodes:
        nodes2kill = vector_of_int()
        for i in range(newnodes.numInternalNodes):
            if mass0[i] == 0.0:
                nodes2kill.append(i)
        if nodes2kill.size() > 0:
            nodes.deleteNodes(nodes2kill)

    # Whew!
    print("Finished resampling nodes: final node count %i." % mpi.allreduce(nodes.numInternalNodes, mpi.SUM))
    return
