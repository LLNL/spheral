#-------------------------------------------------------------------------------
# Dumper class for outputting Spheral++ format data for generating synthetic
# radiographs with Hades.
#-------------------------------------------------------------------------------
import Spheral
import struct
import time

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
# Write a pre-existing Hades format.
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
def hadesDump(integrator,
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
