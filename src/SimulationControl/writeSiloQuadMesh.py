#-------------------------------------------------------------------------------
# Take a set of lattice sampled data (such as from sampleMultipleFields2Lattice)
# and write out to a set of silo files in Quadvar format.
#
# Notes:
#  - Currently only support writing scalar data.  No reason we can't add vectors
#    and tensors down the line.
#-------------------------------------------------------------------------------
def writeSiloQuadMesh(scalar_data,
                      ndim,
                      xmin,
                      xmax,
                      nglobal,
                      filename,
                      dirname = ".",
                      scalar_names = None,
                      materials = [],
                      time = 0.0,
                      cycle = 0,
                      RZ = False):

    import Spheral
    from SpheralCompiledPackages import silo
    import mpi
    import sys, os, struct, bisect
    from operator import mul

    assert ndim in (2,3)

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
            import Spheral2d as sph
        else:
            import Spheral3d as sph

        dx = [(xmax[j] - xmin[j])/nglobal[j] for j in xrange(ndim)]
        ntot = reduce(mul, nglobal)

        # Which dimension should we divide up into?
        jsplit = min(ndim - 1, max(enumerate(nglobal), key = lambda x: x[1])[0])

        # Find the offset to the global lattice numbering on this domain.
        # This is based on knowing the native lattice sampling method stripes the original data
        # accoriding to (i + j*nx + k*nx*ny), and simply divides that 1D serialization sequentially
        # between processors.
        offset = 0
        for sendproc in xrange(mpi.procs):
            n = mpi.bcast(len(vals), root=sendproc)
            if sendproc < mpi.rank:
                offset += n
        if mpi.rank == mpi.procs - 1:
            assert offset + len(vals) == ntot

        # A function to turn an index into the integer lattice coordinates
        def latticeCoords(iglobal):
            return (iglobal % nglobal[0],
                    (iglobal % (nglobal[0]*nglobal[1])) // nglobal[0],
                    iglobal // (nglobal[0]*nglobal[1]))

        # A function to tell us which block to assign a global index to
        slabsperblock = max(1, nglobal[jsplit] // mpi.procs)
        remainder = max(0, nglobal[jsplit] - mpi.procs*slabsperblock)
        islabdomain = [min(nglobal[jsplit], iproc*slabsperblock + min(iproc, remainder)) for iproc in xrange(mpi.procs + 1)]
        #sys.stderr.write("Domain splitting: %s %i %s\n" % (nglobal, jsplit, islabdomain))
        #sys.stderr.write("islabdomain : %s\n" % str(islabdomain))
        def targetBlock(index):
            icoords = latticeCoords(offset + index)
            return bisect.bisect_right(islabdomain, icoords[jsplit]) - 1

        # Build a list of (global_index, value, target_proc) for each of the lattice values.
        id_val_procs = [(offset + i, val, targetBlock(i)) for i, val in enumerate(vals)]
        #sys.stderr.write("id_val_procs : %s\n" % str(id_val_procs))
        #sys.stderr.write("map index -> slab : %s\n" % str([(offset + i, latticeCoords(offset + i), targetBlock(i)) for i in xrange(len(vals))]))
        #sys.stderr.write("id_val_procs : %s\n" % str([(i, tb, latticeCoords(i)) for (i, val, tb) in id_val_procs if i % 100 < 10 and tb != 0]))

        # Send our values to other domains.
        sendreqs, sendvals = [], []
        for iproc in xrange(mpi.procs):
            if iproc != mpi.rank:
                sendvals.append([(i, val) for (i, val, proc) in id_val_procs if proc == iproc])
                sendreqs.append(mpi.isend(sendvals[-1], dest=iproc, tag=100))

        # Now we can build the dang result.
        xminblock, xmaxblock = sph.Vector(*xmin), sph.Vector(*xmax)
        xminblock[jsplit] = xmin[jsplit] + islabdomain[mpi.rank]    *dx[jsplit]
        xmaxblock[jsplit] = xmin[jsplit] + islabdomain[mpi.rank + 1]*dx[jsplit]
        nblock = list(nglobal)
        nblock[jsplit] = islabdomain[mpi.rank + 1] - islabdomain[mpi.rank]
        #sys.stderr.write("nblock : %s\n" % str(nblock))
        newvals = []
        for iproc in xrange(mpi.procs):
            if iproc == mpi.rank:
                recvvals = [(i, val) for (i, val, proc) in id_val_procs if proc == mpi.rank]
            else:
                recvvals = mpi.recv(source=iproc, tag=100)[0]
            newvals += recvvals
        newvals.sort()
        valsblock = sph.vector_of_double()
        for i, val in newvals:
            valsblock.append(val)
        #sys.stderr.write("len(valsblock) = %s\n" % len(valsblock))
        assert len(valsblock) == reduce(mul, nblock)

        # Wait 'til all communication is done.
        for req in sendreqs:
            req.wait()

        # That should be it.
        return valsblock, xminblock, xmaxblock, nblock, jsplit

    #-------------------------------------------------------------------------------
    # Write the master file.
    #-------------------------------------------------------------------------------
    def writeMasterSiloFile(ndim, nblock, jsplit,
                            baseDirectory, baseName, procDir, materials,
                            vars, label, time, cycle):

        nullOpts = silo.DBoptlist()

        # Decide which domains have information.
        if len(vars[0][0]) > 0:
            myvote = mpi.rank + 1
        else:
            myvote = 0
        maxproc = mpi.allreduce(myvote, mpi.MAX)
        assert maxproc <= mpi.procs

        # Pattern for constructing per domain variables.
        domainNamePatterns = [os.path.join(procDir, "domain%i.silo:%%s" % i) for i in xrange(maxproc)]

        # We need each domains nblock info.
        nblocks = [nblock]
        for sendproc in xrange(1, maxproc):
            if mpi.rank == sendproc:
                mpi.send(nblock, dest=0, tag=50)
            if mpi.rank == 0:
                nblocks.append(mpi.recv(source=sendproc, tag=50)[0])

        # Create the master file.
        if mpi.rank == 0:
            fileName = os.path.join(baseDirectory, baseName + ".silo")
            f = silo.DBCreate(fileName, 
                              silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)
            nullOpts = silo.DBoptlist()

            # Write the domain file names and types.
            domainNames = Spheral.vector_of_string([p % "hblk0/hydro_mesh" for p in domainNamePatterns])
            meshTypes = Spheral.vector_of_int([silo.DB_QUADMESH]*maxproc)
            optlist = silo.DBoptlist(1024)
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            assert silo.DBPutMultimesh(f, "hydro_mesh", domainNames, meshTypes, optlist) == 0

            # Write material names.
            if materials:
                material_names = Spheral.vector_of_string([p % "/hblk0/Materials" for p in domainNamePatterns])
                matnames = Spheral.vector_of_string(["void"] + [x.name for x in materials])
                matnos = Spheral.vector_of_int(range(len(materials) + 1))
                assert len(material_names) == maxproc
                assert len(matnames) == len(materials) + 1
                assert len(matnos) == len(materials) + 1
                optlist = silo.DBoptlist(1024)
                assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
                assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
                assert optlist.addOption(silo.DBOPT_MMESH_NAME, "hydro_mesh") == 0
                assert optlist.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, matnames) == 0
                assert optlist.addOption(silo.DBOPT_MATNOS, silo.DBOPT_NMATNOS, matnos) == 0
                assert silo.DBPutMultimat(f, "Materials", material_names, optlist) == 0

            # Write the variables descriptors.
            types = Spheral.vector_of_int([silo.DB_QUADVAR]*maxproc)
            for var, varname in vars:
                domainVarNames = Spheral.vector_of_string()
                for iproc, p in enumerate(domainNamePatterns):
                    domainVarNames.append(p % ("/hblk0/" + varname))
                assert len(domainVarNames) == maxproc
                optlistMV = silo.DBoptlist()
                assert optlistMV.addOption(silo.DBOPT_CYCLE, cycle) == 0
                assert optlistMV.addOption(silo.DBOPT_DTIME, time) == 0
                #assert optlistMV.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_SCALAR) == 0
                assert optlistMV.addOption(silo.DBOPT_BLOCKORIGIN, 0) == 0
                assert optlistMV.addOption(silo.DBOPT_MMESH_NAME, "hydro_mesh") == 0
                assert silo.DBPutMultivar(f, varname, domainVarNames, types, optlistMV) == 0

            # Write the dummy variable "akap_0" to fool Hades into thinking we're actually Hydra or something.
            assert silo.DBPutQuadvar1(f, "akap_0", "hydro_mesh",
                                      Spheral.vector_of_double([0.0]*(ndim*ndim)), Spheral.vector_of_double(),
                                      silo.DB_ZONECENT, Spheral.vector_of_int([ndim]*ndim), nullOpts) == 0

            # Write domain and mesh size info.
            assert silo.DBMkDir(f, "Decomposition") == 0
            assert silo.DBWrite(f, "Decomposition/NumDomains", maxproc) == 0
            assert silo.DBWrite(f, "Decomposition/NumLocalDomains", maxproc) == 0
            assert silo.DBWrite(f, "Decomposition/NumBlocks", 1) == 0
            #assert silo.DBWrite(f, "Decomposition/LocalName", "hblk") == 0
            localDomains = Spheral.vector_of_int(range(maxproc))
            domainFiles = Spheral.vector_of_vector_of_int([Spheral.vector_of_int(range(maxproc))])
            assert silo.DBWrite(f, "Decomposition/LocalDomains", localDomains) == 0
            assert silo.DBWrite(f, "DomainFiles", domainFiles) == 0

            for iproc in xrange(maxproc):
                assert silo.DBMkDir(f, "Decomposition/gmap%i" % iproc) == 0
                stuff = Spheral.vector_of_int([0]*12)
                for jdim in xrange(ndim):
                    stuff[6+jdim] = nblocks[iproc][jdim]
                if iproc in (0, maxproc-1):
                    assert silo.DBWrite(f, "Decomposition/gmap%i/NumNeighbors" % iproc, 1) == 0
                else:
                    assert silo.DBWrite(f, "Decomposition/gmap%i/NumNeighbors" % iproc, 2) == 0
                assert silo.DBWrite(f, "Decomposition/gmap%i/gmap" % iproc, stuff) == 0

        # Close the file.
        if mpi.rank == 0:
            assert silo.DBClose(f) == 0
            del f

        return maxproc

    #-------------------------------------------------------------------------------
    # Write the domain file.
    #-------------------------------------------------------------------------------
    def writeDomainSiloFile(ndim, maxproc,
                            baseDirectory, baseName, procDir,
                            materials, vars,
                            xminblock, xmaxblock, nblock, jsplit,
                            label, time, cycle,
                            pretendRZ):

        assert jsplit < ndim

        # Make sure the directories are there.
        if mpi.rank == 0:
            for iproc in xrange(maxproc):
                pth = os.path.join(baseDirectory, procDir)
                if not os.path.exists(pth):
                    os.makedirs(pth)
        mpi.barrier()

        # Is there anything to do?
        if mpi.rank < maxproc:
            numZones = 1
            numNodes = 1
            nblocknodes = list(nblock)
            for i, x in enumerate(nblock):
                numZones *= x
                numNodes *= x + 1
                nblocknodes[i] = x + 1
            assert numZones > 0

            # Make a vector<int> version of nblock
            nblock_vec = Spheral.vector_of_int(nblock)

            # Create the file.
            fileName = os.path.join(baseDirectory, procDir, "domain%i.silo" % mpi.rank)
            f = silo.DBCreate(fileName, 
                              silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)
            nullOpts = silo.DBoptlist()

            # Make the hblk0 directory.
            assert silo.DBMkDir(f, "hblk0") == 0

            # Write the domain mesh.
            coords = Spheral.vector_of_vector_of_double([Spheral.vector_of_double()]*ndim)
            for jdim in xrange(ndim):
                coords[jdim] = Spheral.vector_of_double([0.0]*nblocknodes[jdim])
                dx = (xmaxblock[jdim] - xminblock[jdim])/nblock[jdim]
                for i in xrange(nblocknodes[jdim]):
                    coords[jdim][i] = xminblock[jdim] + i*dx
            optlist = silo.DBoptlist()
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            if pretendRZ:
                assert optlist.addOption(silo.DBOPT_COORDSYS, silo.DB_CYLINDRICAL) == 0
            else:
                assert optlist.addOption(silo.DBOPT_COORDSYS, silo.DB_CARTESIAN) == 0
            assert silo.DBPutQuadmesh(f, "hblk0/hydro_mesh", coords, optlist) == 0

            # Write materials.
            if materials:
                matnos = Spheral.vector_of_int(range(len(materials)+1))
                assert len(matnos) == len(materials) + 1
                matlist = Spheral.vector_of_int([0]*numZones)
                matnames = Spheral.vector_of_string(["void"])
                for imat, nodeList in enumerate(materials):
                    for i in xrange(numZones):
                        if vars[0][0][i] > 0.0:
                            matlist[i] = imat + 1
                    matnames.append(nodeList.name)
                assert len(matlist) == numZones
                assert len(matnames) == len(materials) + 1
                matOpts = silo.DBoptlist(1024)
                assert matOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
                assert matOpts.addOption(silo.DBOPT_DTIME, time) == 0
                assert matOpts.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, matnames) == 0
                assert silo.DBPutMaterial(f, "hblk0/Materials", "hydro_mesh", matnos, matlist, nblock_vec,
                                          Spheral.vector_of_int(), Spheral.vector_of_int(), Spheral.vector_of_int(), Spheral.vector_of_double(),
                                          matOpts) == 0

            # Write the field variables.
            varOpts = silo.DBoptlist(1024)
            assert varOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert varOpts.addOption(silo.DBOPT_DTIME, time) == 0
            for var, varname in vars:
                assert len(var) == numZones
                assert silo.DBPutQuadvar1(f, "hblk0/" + varname, "hydro_mesh", var, 
                                          Spheral.vector_of_double(), silo.DB_ZONECENT, nblock_vec, varOpts) == 0

            # That's it.
            assert silo.DBClose(f) == 0
            del f

        return

    #---------------------------------------------------------------------------
    # And finally the main routine...
    #---------------------------------------------------------------------------

    # Output field names for file.
    if scalar_names:
        assert len(scalar_names) == len(scalar_data)
    else:
        scalar_names = ["scalar%i" % i for i in xrange(len(scalar_names))]

    # Shuffle the scalar data into the block array structure required by silo.
    scalar_blocks = []
    for var in scalar_data:
        varsamp, xminblock, xmaxblock, nblock, jsplit = shuffleIntoBlocks(ndim, var, xmin, xmax, nglobal)
        scalar_blocks.append(varsamp)

    procDir = "blocks_cycle%i" % cycle

    # Write the files.
    maxproc = writeMasterSiloFile(ndim = ndim,
                                  nblock = nblock,
                                  jsplit = jsplit,
                                  baseDirectory = dirname,
                                  baseName = filename,
                                  procDir = procDir,
                                  materials = materials,
                                  vars = zip(scalar_blocks, scalar_names),
                                  label = "Spheral++ cartesian sampled output",
                                  time = time,
                                  cycle = cycle)
    writeDomainSiloFile(ndim = ndim,
                        jsplit = jsplit,
                        maxproc = maxproc,
                        baseDirectory = dirname,
                        baseName = filename,
                        procDir = procDir,
                        materials = materials,
                        vars = zip(scalar_blocks, scalar_names),
                        xminblock = xminblock,
                        xmaxblock = xmaxblock,
                        nblock = nblock,
                        label = "Spheral++ cartesian sampled output",
                        time = time,
                        cycle = cycle,
                        pretendRZ = RZ)

