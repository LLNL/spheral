#-------------------------------------------------------------------------------
# Dumper class for outputting Spheral++ format data for generating synthetic
# radiographs with Hades.
#-------------------------------------------------------------------------------
import Spheral
from SpheralModules import silo
from SpheralModules.silo import SiloAttributes as SA
import mpi
import struct, time, os, bisect

#-------------------------------------------------------------------------------
# Write a silo file resampling to a fixed cartesian mesh for the density.
#-------------------------------------------------------------------------------
def hadesDump(integrator,
              nsample,
              xmin,
              xmax,
              W,
              baseFileName,
              baseDirectory = ".",
              procDirBaseName = "proc-%06i",
              mask = None,
              nodeLists = None):

    # Currently suppport 2D and 3D.
    db = integrator.dataBase()
    if db.nDim == 2:
        import Spheral2d as sph
    elif db.nDim == 3:
        import Spheral3d as sph
    else:
        raise RuntimeError, "hadesDump ERROR: must be 2D or 3D"

    # Prepare to time how long this takes.
    t0 = time.clock()

    # Get the set of material names we're going to write.
    if materials is None:
        materials = list(db.fluidNodeLists())

    # HACK!  We are currently restricting to writing single material output!
    assert len(materials) == 1

    # Make sure the output directory exists.
    if mpi.rank == 0 and not os.path.exists(baseDirectory):
        try:
            os.makedirs(baseDirectory)
        except:
            raise RuntimeError, "Cannot create output directory %s" % baseDirectory
    mpi.barrier()

    # Sample the density.
    ntot = nsample[0]
    for j in xrange(2, db.nDim):
        ntot *= nsample[j]
    for nodes in materials:
        print "hadesDump: sampling density for %s..." % nodes.name
        r = Spheral.VectorFieldList()
        H = Spheral.SymTensorFieldList()
        rho = Spheral.ScalarFieldList()
        r.appendField(nodes.positions())
        w.appendField(nodes.weight())
        H.appendField(nodes.Hfield())
        rho.appendField(nodes.massDensity())

        w = Spheral.ScalarFieldList()
        w.copyFields()
        w.appendField(Spheral.ScalarField("weight", nodes, 1.0))

        fieldListSet = Spheral.FieldListSet()
        fieldListSet.ScalarFieldLists.append(rho)
        localMask = Spheral.IntFieldList()
        if mask is None:
            localMask.copyFields()
            localMask.appendField(Spheral.IntField("mask", nodes, 1))
        else:
            localMask.appendField(mask.fieldForNodeList(nodes))

        scalar_samples = Spheral.vector_of_vector_of_double()
        vector_samples = Spheral.vector_of_vector_of_Vector()
        tensor_samples = Spheral.vector_of_vector_of_Tensor()
        symTensor_samples = Spheral.vector_of_vector_of_SymTensor()
        nsample_vec = Spheral.vector_of_int(db.nDim)
        for i in xrange(db.nDim):
            nsample_vec[i] = nsample[i]

        Spheral.sampleMultipleFields2Lattice(fieldListSet,
                                             r, w, H, localMask,
                                             W,
                                             xmin, xmax,
                                             nsample_vec,
                                             scalar_samples,
                                             vector_samples,
                                             tensor_samples,
                                             symTensor_samples)
        print "Generated %i scalar fields" % len(scalar_samples)

        # Rearrange the sampled data into rectangular blocks due to Silo's quad mesh limitations.
        rhosamp, xminblock, xmaxblock, nblock = shuffleIntoBlocks(ndim, scalar_fields[0], xmin, xmax, nsample)
        assert mpi.allreduce(len(rhosamp), mpi.SUM) == ntot

    # Write the master file.
    writeMasterSiloFile(ndim = db.nDim,
                        baseDirectory = baseDirectory,
                        baseName = baseFileName,
                        procDirBaseName = procDirBaseName,
                        nodeLists = materials,
                        rhosamp = rhosamp,
                        label = "Spheral++ cartesian sampled output",
                        time = integrator.currentTime,
                        cycle = integrator.currentCycle)

    # Write the process files.
    writeDomainSiloFile(ndim = db.nDim,
                        dirName = proc
                        )

    mpi.barrier()
    print "hadesDump finished: required %0.2f seconds" % (time.clock() - t0)
    return

#-------------------------------------------------------------------------------
# Rearrange the distributed data into rectangular blocks for each domain due
# to Silo's restrictions for quad meshes.
# Returns:
#   valsblock : array of local values arranged in a block
#   xminblock : min coordinates of block
#   xmaxblock : max coordinates of block
#   nblock    : ndim dimensioned array, length of sampled block in each dimension
#
# This method simplifies by always cutting into slabs.  This may result in
# some domains having zero values though, if there are more processors than
# lattice slabs in the chosen direction.
#-------------------------------------------------------------------------------
def shuffleIntoBlocks(ndim, vals, xmin, xmax, nglobal):

    if ndim == 2:
        from Spheral2d import *
    else:
        from Spheral3d import *

    # Which dimension are cutting up into?
    jmax = max(enumerate(nglobal), key = lambda x: x[1])[0]
    dx = [(xmax[j] - xmin[j])/nglobal[j] for j in xrange(ndim)]
    ntot = 1
    for x in nglobal:
        ntot *= x

    # Find the offset to the global lattice numbering on this domain
    offset = 0
    for sendproc in xrange(mpi.procs):
        n = mpi.bcast(len(vals), root=sendproc)
        if sendproc < mpi.rank:
            offset += n
    if mpi.rank == mpi.procs - 1:
        assert offset + len(vals) == ntot

    # A function to tell us which block to assign a global index to
    nperslab = 1
    for j in xrange(ndim):
        if j != jmax:
            nperslab *= nglobal[j]
    slabsperblock = max(1, nglobal[jmax] // mpi.procs)
    remainder = max(0, nglobal[jmax] - mpi.procs*slabsperblock)
    islabdomain = [min(nglobal[jmax], iproc*slabsperblock + min(iproc, remainder)) for iproc in xrange(mpi.procs)]
    def targetBlock(index):
        index = offset + index
        if ndim == 2:
            id = (index % nglobal[0],
                  index // nglobal[0],
                  0)
        else:
            id = ((index % (nglobal[0]*nglobal[1])) % nglobal[0],
                  (index % (nglobal[0]*nglobal[1])) // nglobal[0],
                  index // (nglobal[0]*nglobal[1]))
        return bisect.bisect(islabdomain, id[jmax])

    # Now we can build the dang result.
    xminblock, xmaxblock = Vector(xmin), Vector(xmax)
    xminblock[jmax] = islabdomain[mpi.rank]    *dx[jmax]
    xmaxblock[jmax] = islabdomain[mpi.rank + 1]*dx[jmax]
    nblock = list(nglobal)
    nblock[jmax] = islabdomain[mpi.rank + 1] - islabdomain[mpi.rank]
    valsblock = []
    for i, vali in enumerate(vals):
        

#-------------------------------------------------------------------------------
# Write the master file.
#-------------------------------------------------------------------------------
def writeMasterSiloFile(ndim, baseDirectory, baseName, procDirBaseName, nodeLists,
                        rhosamp, label, time, cycle):

    nullOpts = silo.DBoptlist()

    # Pattern for constructing per domain variables.
    domainNamePatterns = [os.path.join(procDirBaseName % i, baseName + ".silo:%s") for i in xrange(mpi.procs)]

    # Create the master file.
    if mpi.rank == 0:
        fileName = os.path.join(baseDirectory, baseName + ".silo")
        f = silo.DBCreate(fileName, 
                          SA._DB_CLOBBER, SA._DB_LOCAL, label, SA._DB_HDF5)

        # Write the domain file names and types.
        domainNames = vector_of_string()
        meshTypes = vector_of_int(mpi.procs, SA._DB_QUAD_RECT)
        for p in domainNamePatterns:
            domainNames.append(p % "MESH")
        optlist = silo.DBoptlist(1024)
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
        assert silo.DBPutMultimesh(f, "MMESH", domainNames, meshTypes, optlist) == 0

        # Write material names.
        material_names = vector_of_string()
        matnames = vector_of_string()
        matnos = vector_of_int()
        for p in domainNamePatterns:
            material_names.append(p % "material")
        for i, name in enumerate([x.name for x in nodeLists]):
            matnames.append(name)
            matnos.append(i)
        assert len(material_names) == mpi.procs
        assert len(matnames) == len(nodeLists)
        assert len(matnos) == len(nodeLists)
        optlist = silo.DBoptlist(1024)
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
        assert optlist.addOption(SA._DBOPT_MATNAMES, SA._DBOPT_NMATNOS, matnames) == 0
        assert optlist.addOption(SA._DBOPT_MATNOS, SA._DBOPT_NMATNOS, matnos) == 0
        assert silo.DBPutMultimat(f, "MMATERIAL", material_names, optlist) == 0
        
        # Write the variables descriptors.
        # We currently hardwire for the single density variable.
        name = "mass_density"
        types = vector_of_int(mpi.procs, SA._DB_QUADVAR)
        domainVarNames = vector_of_string()
        nlocalvals = len(rhosamp)
        for iproc, p in enumerate(domainNamePatterns):
            nvals = mpi.bcast(nlocalvals, root=iproc)
            if nvals > 0:
                domainVarNames.append(p % name)
            else:
                domainVarNames.append("EMPTY")
        assert len(domainVarNames) == mpi.procs
        optlistMV = silo.DBoptlist()
        assert optlistMV.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlistMV.addOption(SA._DBOPT_DTIME, time) == 0
        assert optlistMV.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_SCALAR) == 0
        if mpi.rank == 0:
            assert silo.DBPutMultivar(f, name, domainVarNames, types, optlistMV) == 0

    # That's it.
    if mpi.rank == 0:
        assert silo.DBClose(f) == 0
        del f

    return

#-------------------------------------------------------------------------------
# Write the domain file.
#-------------------------------------------------------------------------------
def writeDomainSiloFile(ndim, procDirBaseName, mesh, index2zone, label, nodeLists, time, cycle, fieldwad,
                        pretendRZ, nodeArrays, zoneArrays, faceArrays,
                        meshType = SA._DB_UCDMESH):

    # Is there anything to do?
    numZones = len(mesh.cells)
    if numZones > 0:

        # Create the file.
        fileName = os.path.join(dirName, "domain%i.silo" % mpi.rank)
        db = silo.DBCreate(fileName, 
                           SA._DB_CLOBBER, SA._DB_LOCAL, label, SA._DB_HDF5)
        nullOpts = silo.DBoptlist()

        # Determine our dimensionality
        if isinstance(mesh, polytope.Tessellation2d):
            nDim = 2
        else:
            assert isinstance(mesh, polytope.Tessellation3d)
            nDim = 3

        # Write a Polygonal zone list.
        zonelistName = { 2 : "zonelist",
                         3 : "PHzonelist" }

        if nDim == 2:
        
            # Read out the zone nodes.  We rely on these already being arranged
            # counter-clockwise.
            zoneNodes = vector_of_vector_of_int()
            shapesize = vector_of_int()
            for zoneID in xrange(numZones):
                zone = mesh.cells[zoneID]
                nodes = vector_of_int()
                for iface in zone:
                    if iface < 0:
                        nodes.append(mesh.faces[~iface][1])
                    else:
                        nodes.append(mesh.faces[iface][0])
                zoneNodes.append(nodes)
                shapesize.append(len(nodes))
            assert len(zoneNodes) == numZones
            assert len(shapesize) == numZones
        
            assert silo.DBPutZonelist2(db, zonelistName[nDim], nDim, zoneNodes, 0, 0,
                                       vector_of_int(numZones, SA._DB_ZONETYPE_POLYGON),
                                       shapesize,
                                       vector_of_int(numZones, 1),
                                       nullOpts) == 0
        
        # Write a Polyhedral zone list.
        if nDim == 3:
        
            # Construct the face-node lists.
            numFaces = len(mesh.faces)
            faceNodes = vector_of_vector_of_int(numFaces)
            for iface in xrange(numFaces):
                for j in xrange(len(mesh.faces[iface])):
                    faceNodes[iface].append(mesh.faces[iface][j])
                assert len(faceNodes[iface]) == len(mesh.faces[iface])
            assert len(faceNodes) == numFaces
        
            # Construct the zone-face list.  We use the ones complement of a face ID
            # to indicate that face needs to be reversed in reference to this zone.
            # This is the same convention as polytope, so just copy it.
            zoneFaces = mesh.cells
            assert len(zoneFaces) == numZones
        
            assert silo.DBPutPHZonelist(db, zonelistName[nDim], faceNodes, zoneFaces, 0, (numZones - 1), nullOpts) == 0
        
        # Construct the mesh node coordinates.
        assert len(mesh.nodes) % nDim == 0
        numNodes = len(mesh.nodes)/nDim
        coords = vector_of_vector_of_double(nDim, vector_of_double(numNodes))
        for nodeID in xrange(numNodes):
            for idim in xrange(nDim):
                coords[idim][nodeID] = mesh.nodes[nDim*nodeID + idim]
        assert len(coords) == nDim
        
        # Write the mesh itself.
        meshOpts = silo.DBoptlist(1024)
        assert meshOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert meshOpts.addOption(SA._DBOPT_DTIME, time) == 0
        assert meshOpts.addOption(SA._DBOPT_COORDSYS, SA._DB_CARTESIAN) == 0
        assert meshOpts.addOption(SA._DBOPT_NSPACE, nDim) == 0
        assert meshOpts.addOption(SA._DBOPT_TV_CONNECTIVITY, 1) == 0
        if nDim == 2:
            if pretendRZ:
                assert meshOpts.addOption(SA._DBOPT_COORDSYS, SA._DB_CYLINDRICAL) == 0
                assert meshOpts.addOption(SA._DBOPT_XLABEL, "z") == 0
                assert meshOpts.addOption(SA._DBOPT_YLABEL, "r") == 0
            assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, zonelistName[nDim], "NULL", meshOpts) == 0
        else:
            assert meshOpts.addOption(SA._DBOPT_PHZONELIST, zonelistName[nDim]) == 0
            assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, "NULL", "NULL", meshOpts) == 0
        
        # Write materials.
        if nodeLists:
            matnos = vector_of_int()
            for i in xrange(len(nodeLists)):
                matnos.append(i)
            assert len(matnos) == len(nodeLists)
            matlist = vector_of_int(numZones)
            matnames = vector_of_string()
            offset = 0
            for (nodeList, imat) in zip(nodeLists, xrange(len(nodeLists))):
                for i in xrange(nodeList.numInternalNodes):
                    for j in index2zone[offset + i]:
                        matlist[j] = imat
                matnames.append(nodeList.name)
                offset += nodeList.numInternalNodes
            assert len(matlist) == numZones
            assert len(matnames) == len(nodeLists)
            matOpts = silo.DBoptlist(1024)
            assert matOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
            assert matOpts.addOption(SA._DBOPT_DTIME, time) == 0
            assert matOpts.addOption(SA._DBOPT_MATNAMES, SA._DBOPT_NMATNOS, matnames) == 0
            assert silo.DBPutMaterial(db, "MATERIAL", "MESH", matnos, matlist,
                                      vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                                      matOpts) == 0
            assert silo.DBPutMaterial(db, "PointMATERIAL", "PointMESH", matnos, matlist,
                                      vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                                      matOpts) == 0
        
            # Write the variable descriptions for non-scalar variables (vector and tensors).
            writeDefvars(db, fieldwad)
        
            # Write the field components.
            centering = SA._DB_ZONECENT
            varOpts = silo.DBoptlist(1024)
            assert varOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
            assert varOpts.addOption(SA._DBOPT_DTIME, time) == 0
            for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
                for subname, vals in subvars:
                    if len(vals) > 0:
                        if isinstance(vals, vector_of_double):
                            assert silo.DBPutUcdvar1(db, subname, "MESH", vals, vector_of_double(), centering, varOpts) == 0
                        elif isinstance(vals, vector_of_int):
                            assert silo.DBPutUcdvar1(db, subname, "MESH", vals, vector_of_int(), centering, varOpts) == 0
        
        # Write the set of neighbor domains.
        thpt = vector_of_vector_of_int()
        thpt.append(vector_of_int(1, len(mesh.neighborDomains)))
        thpt.append(vector_of_int())
        for i in mesh.neighborDomains:
            thpt[-1].append(i)
        elemNames = vector_of_string()
        elemNames.append("num neighbor domains")
        elemNames.append("neighbor domains")
        assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR_NUMS", elemNames, thpt, nullOpts) == 0
        
        # Write the shared nodes for each neighbor domain.
        sharedNodes = mesh.sharedNodes
        for ineighborDomain in xrange(len(mesh.neighborDomains)):
            nodes = vector_of_int()
            for i in sharedNodes[ineighborDomain]:
                nodes.append(i)
            assert len(nodes) == len(sharedNodes[ineighborDomain])
            assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR%i" % ineighborDomain,
                                           vector_of_string(1, "shared_nodes"),
                                           vector_of_vector_of_int(1, nodes),
                                           nullOpts) == 0
        
        # If requested, write out annotations for the nodes, zones, and faces.
        if (not (nodeArrays is None) or
            not (zoneArrays is None) or
            not (faceArrays is None)):
            names = vector_of_string()
            values = vector_of_vector_of_int()
            if not (nodeArrays is None):
                for pair in nodeArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_node")
                        values.append(vector_of_int())
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            if not (zoneArrays is None):
                for pair in zoneArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_zone")
                        values.append(vector_of_int())
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            if not (faceArrays is None):
                for pair in faceArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_face")
                        values.append(vector_of_int())
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            assert len(names) == len(values)
            if len(names) > 0:
                assert silo.DBPutCompoundarray(db, "ANNOTATION_INT", names, values, nullOpts) == 0
        
        # Write the point mesh.
        if nodeLists:
            ntot = sum([n.numInternalNodes for n in nodeLists])
            coords = vector_of_vector_of_double(nDim)
            for nodes in nodeLists:
                pos = nodes.positions().internalValues()
                n = len(pos)
                for j in xrange(nDim):
                    for i in xrange(n):
                        coords[j].append(pos[i][j])
            for j in xrange(nDim):
                assert len(coords[j]) == ntot
         
            # Write the Pointmesh.
            meshOpts = silo.DBoptlist(1024)
            assert meshOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
            assert meshOpts.addOption(SA._DBOPT_DTIME, time) == 0
            assert silo.DBPutPointmesh(db, "PointMESH", coords, meshOpts) == 0

        # That's it.
        assert silo.DBClose(db) == 0
        del db

    return

#-------------------------------------------------------------------------------
# Our original proposition for a Hades file format.
#-------------------------------------------------------------------------------
## hadesHeader = """
## # Spheral++ -> Hades file
## # Format:
## #   N (Num materials)
## #   xmin ymin zmin
## #   xmax ymax zmax
## #   nx ny nz
## #   (num isotopes 1) ZA frac ZA frac ...
## #   ...
## #   (num isotopes N) ZA frac ZA frac ...
## #   mass density 1
## #   ...
## #   mass density N
## # Mass densities are written as index1 rho1 index2 rho2 ...
## """
def hadesDump0(integrator,
               nsample,
               xmin,
               xmax,
               W,
               isotopes,
               baseFileName,
               baseDirectory = ".",
               dumpGhosts = False,
               materials = "all"):

    # We currently only support 3-D.
    assert isinstance(integrator, Spheral.Integrator3d)
    assert len(nsample) == 3
    assert isinstance(xmin, Spheral.Vector3d)
    assert isinstance(xmax, Spheral.Vector3d)
    assert isinstance(W, Spheral.TableKernel3d)
    for x in isotopes:
        for xx in x:
            assert len(xx) == 2

    # Prepare to time how long this takes.
    t0 = time.clock()

    # Extract the data base.
    db = integrator.dataBase()

    # If requested, set ghost node info.
    if dumpGhosts and not integrator is None:
        state = Spheral.State3d(db, integrator.physicsPackages())
        derivs = Spheral.StateDerivatives3d(db, integrator.physicsPackages())
        integrator.setGhostNodes()
        integrator.applyGhostBoundaries(state, derivs)

    # Get the set of material names we're going to write.
    if materials == "all":
        materials = [n for n in db.fluidNodeLists()]
    assert len(materials) == len(isotopes)

    # Make sure the output directory exists.
    import mpi
    import os
    if mpi.rank == 0 and not os.path.exists(baseDirectory):
        try:
            os.makedirs(baseDirectory)
        except:
            raise "Cannot create output directory %s" % baseDirectory
    mpi.barrier()

    # Open a file for the output.
    currentTime = integrator.currentTime
    currentCycle = integrator.currentCycle
    filename = baseDirectory + "/" + baseFileName + "-time=%g-cycle=%i.hades" % (currentTime, currentCycle)

    if mpi.rank == 0:
        f = open(filename, "wb")

        # Write the header info.
        #f.write(hadesHeader)
        f.write(struct.pack("I", len(materials)))
        f.write(struct.pack("ddd", *tuple(xmin.elements())))
        f.write(struct.pack("ddd", *tuple(xmax.elements())))
        f.write(struct.pack("III", *nsample))
        for materialIsotopes in isotopes:
            f.write(struct.pack("I", len(materialIsotopes)))
            for iso in materialIsotopes:
                f.write(struct.pack("Id", *iso))

    # For each material, sample the mass density and write it out.
    ntot = nsample[0]*nsample[1]*nsample[2]
    for nodes in materials:
        r = Spheral.VectorFieldList3d()
        w = Spheral.ScalarFieldList3d()
        H = Spheral.SymTensorFieldList3d()
        rho = Spheral.ScalarFieldList3d()
        r.appendField(nodes.positions())
        w.appendField(nodes.weight())
        H.appendField(nodes.Hfield())
        rho.appendField(nodes.massDensity())
        fieldListSet = Spheral.FieldListSet3d()
        fieldListSet.ScalarFieldLists.append(rho)
        rhosamp = Spheral.sampleMultipleFields2LatticeMash(fieldListSet,
                                                           r, w, H,
                                                           W,
                                                           xmin, xmax,
                                                           nsample)[0][0][1]
        assert mpi.allreduce(len(rhosamp), mpi.SUM) == ntot
        icum = 0
        for sendProc in xrange(mpi.procs):
            valsproc = [(i, x) for (i, x) in zip(range(ntot), rhosamp) if x > 0.0]
            vals = mpi.bcast(valsproc, sendProc)
            if mpi.rank == 0:
                f.write(struct.pack("I", len(vals)))
                for i, x in vals:
                    f.write(struct.pack("id", i + icum, x))
            icum += len(vals)

    if mpi.rank == 0:
        # Close the file and we're done.
        f.close()

    mpi.barrier()
    print "hadesDump finished: required %0.2f seconds" % (time.clock() - t0)

    return

#-------------------------------------------------------------------------------
# Write a pre-existing Hades format. (Another defunct option)
#
# This format writes four files:
# 
# <basename>.spr (ASCII)
#   3     # (I), number of dimensions
#   nx    # (I), numer in x direction
#   xmin  # (F), minimum x value
#   lx    # (F), cell length in x
#   ny    # (I), numer in y direction
#   ymin  # (F), minimum y value
#   ly    # (F), cell length in y
#   nz    # (I), numer in z direction
#   zmin  # (F), minimum z value
#   lz    # (F), cell length in z
#   3     # (I), HADES flag indicating type
#
# <basename>.sdt (binary)
# <float data, nx x ny x nz>  # Full density data as float*4 array
#
# <basename)>_mat.sdt  (binary)
# <int data, nx x ny x nz>    # Integer flags indicating material in each cell
#
# isos.mat       (ASCII)
#  isofracs(1) = ZA frac ZA frac ....   # isotopics for material 1
#  isofracs(2) = ZA frac ZA frac ....   # isotopics for material 2
#  ...
#  isofracs(n) = ZA frac ZA frac ....   # isotopics for material n
#-------------------------------------------------------------------------------
def hadesDump1(integrator,
              nsample,
              xmin,
              xmax,
              W,
              isotopes,
              baseFileName,
              baseDirectory = ".",
              mask = None,
              dumpGhosts = True,
              materials = "all"):

    # We currently only support 3-D.
    assert isinstance(integrator, Spheral.Integrator3d)
    assert len(nsample) == 3
    assert isinstance(xmin, Spheral.Vector3d)
    assert isinstance(xmax, Spheral.Vector3d)
    assert isinstance(W, Spheral.TableKernel3d)
    for x in isotopes:
        for xx in x:
            assert len(xx) == 2

    # Prepare to time how long this takes.
    t0 = time.clock()

    # Extract the data base.
    db = integrator.dataBase()

    # If requested, set ghost node info.
    if dumpGhosts and not integrator is None:
        state = Spheral.State3d(db, integrator.physicsPackages())
        derivs = Spheral.StateDerivatives3d(db, integrator.physicsPackages())
        integrator.setGhostNodes()
        integrator.applyGhostBoundaries(state, derivs)

    # Get the set of material names we're going to write.
    if materials == "all":
        materials = [n for n in db.fluidNodeLists()]
    assert len(materials) == len(isotopes)

    # HACK!  We are currently restricting to writing single material output!
    assert len(materials) == 1

    # Make sure the output directory exists.
    import mpi
    import os
    if mpi.rank == 0 and not os.path.exists(baseDirectory):
        try:
            os.makedirs(baseDirectory)
        except:
            raise "Cannot create output directory %s" % baseDirectory
    mpi.barrier()

    # Write the density header file.
    print "hadesDump: writing density header..."
    if mpi.rank == 0:
        filename = os.path.join(baseDirectory, baseFileName + ".spr")
        f = open(filename, "w")
        f.write("3\n")
        for i in xrange(3):
            f.write("%i\n" % nsample[i])
            f.write("%f\n" % xmin(i))
            f.write("%f\n" % ((xmax(i) - xmin(i))/nsample[i]))
        f.write("3\n")
        f.close()
    mpi.barrier()

    # Sample the density.
    ntot = nsample[0]*nsample[1]*nsample[2]
    for nodes in materials:
        print "hadesDump: sampling density for %s..." % nodes.name
        r = Spheral.VectorFieldList3d()
        H = Spheral.SymTensorFieldList3d()
        rho = Spheral.ScalarFieldList3d()
        r.appendField(nodes.positions())
        w.appendField(nodes.weight())
        H.appendField(nodes.Hfield())
        rho.appendField(nodes.massDensity())

        w = Spheral.ScalarFieldList3d()
        w.copyFields()
        w.appendField(Spheral.ScalarField3d("weight", nodes, 1.0))

        fieldListSet = Spheral.FieldListSet3d()
        fieldListSet.ScalarFieldLists.append(rho)
        localMask = Spheral.IntFieldList3d()
        if mask is None:
            localMask.copyFields()
            localMask.appendField(Spheral.IntField3d("mask", nodes, 1))
        else:
            localMask.appendField(mask.fieldForNodeList(nodes))

        scalar_samples = Spheral.vector_of_vector_of_double()
        vector_samples = Spheral.vector_of_vector_of_Vector3d()
        tensor_samples = Spheral.vector_of_vector_of_Tensor3d()
        symTensor_samples = Spheral.vector_of_vector_of_SymTensor3d()
        nsample_vec = Spheral.vector_of_int(3)
        for i in xrange(3):
            nsample_vec[i] = nsample[i]

        Spheral.sampleMultipleFields2Lattice3d(fieldListSet,
                                               r, w, H, localMask,
                                               W,
                                               xmin, xmax,
                                               nsample_vec,
                                               scalar_samples,
                                               vector_samples,
                                               tensor_samples,
                                               symTensor_samples)
        print "Generated %i scalar fields" % len(scalar_samples)
        rhosamp = scalar_fields[0]
        nlocal = len(rhosamp)
        assert mpi.allreduce(nlocal, mpi.SUM) == ntot

        print "hadesDump: writing density for %s..." % nodes.name
        filename = os.path.join(baseDirectory, baseFileName + ".sdt")
        for sendProc in xrange(mpi.procs):
            if mpi.rank == sendProc:
                f = open(filename, "ab")
                f.write(struct.pack(nlocal*"f", *tuple(rhosamp)))
                f.close()
            mpi.barrier()

    # Write the material arrays.
    print "hadesDump: writing material flags..."
    filename = os.path.join(baseDirectory, baseFileName + "_mat.sdt")
    for sendProc in xrange(mpi.procs):
        if mpi.rank == sendProc:
            f = open(filename, "ab")
            f.write(struct.pack(nlocal*"i", *(nlocal*(1,))))
            f.close()
        mpi.barrier()

    # Write the isotopes.
    print "hadesDump: writing isotopics..."
    if mpi.rank == 0:
        filename = os.path.join(baseDirectory, "isos.mat")
        f = open(filename, "w")
        i = 0
        for isofracs in isotopes:
            f.write("isofrac(%i) =" % i)
            for (iso, frac) in isofracs:
                f.write(" %i %f" % (iso, frac))
            f.write("\n")
            i += 1
        f.close()
    mpi.barrier()

    mpi.barrier()
    print "hadesDump finished: required %0.2f seconds" % (time.clock() - t0)
    return

