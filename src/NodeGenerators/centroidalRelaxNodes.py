from math import *
import mpi

import Spheral
try:
    from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
except:
    print "centroidalRelaxNoddes unable to import SpheralVoronoiSiloDump -- no tessellation output supported."
    SpheralVoronoiSiloDump = None

#-------------------------------------------------------------------------------
# Centroidally (in mass) relax points allowing a linear density gradient.
#-------------------------------------------------------------------------------
def centroidalRelaxNodes(nodeListsAndBounds,
                         W,
                         rho,
                         gradrho = None,
                         boundaries = [],
                         maxIterations = 100,
                         fracTol = 1.0e-3,
                         correctionOrder = Spheral.LinearOrder,
                         centroidFrac = 0.5,
                         tessellationFileName = None):

    # Did we get passed a function or a constant for the density?
    if type(rho) is float:
        def rhofunc(posi):
            return rho
    else:
        rhofunc = rho

    # What about the gradrho?  Did we get passed anything?
    if gradrho is None:
        gradrhofunc = None
    else:
        if type(gradrho) is float:
            def gradrhofunc(posi):
                return gradrho
        else:
            gradrhofunc = gradrho

    # Decide on our dimensionality and import the appropriate aliases.
    assert (isinstance(W, Spheral.TableKernel1d) or
            isinstance(W, Spheral.TableKernel2d) or
            isinstance(W, Spheral.TableKernel3d))
    if isinstance(W, Spheral.TableKernel1d):
        import Spheral1d as sph
        ndim = 1
    elif isinstance(W, Spheral.TableKernel2d):
        import Spheral2d as sph
        ndim = 2
    else:
        import Spheral3d as sph
        ndim = 3

    # Split out the NodeLists and bounding volumes (if available), depending on what was passed.
    bounds = sph.vector_of_FacetedVolume()
    holes = sph.vector_of_vector_of_FacetedVolume()
    for x in nodeListsAndBounds:
        if type(nodeListsAndBounds[0]) is tuple:
            bounds.resize(len(nodeListsAndBounds))
            holes.resize(len(nodeListsAndBounds))
    if len(bounds) > 0:
        nodeLists = []
        for i, xtup in enumerate(nodeListsAndBounds):
            if type(xtup) is tuple:
                nodeLists.append(xtup[0])
                assert len(xtup) in (2,3)
                bounds[i] = xtup[1]
                if len(xtup) == 3:    # Check for holes
                    assert type(xtup[2]) is list
                    for x in xtup[2]:
                        holes[i].append(x)
            else:
                nodeLists.append(xtup)
    else:
        nodeLists = nodeListsAndBounds

    # Build a local DataBase.
    db = sph.DataBase()
    for nodes in nodeLists:
        db.appendNodeList(nodes)

    # Get references to state in the NodeLists.
    pos = db.fluidPosition
    H = db.fluidHfield
    mass = db.fluidMass
    rhof = db.fluidMassDensity

    # Prepare the storage for the point-wise fields.
    gradRhof = db.newFluidVectorFieldList(sph.Vector.zero, "mass density gradient")
    surfacePoint = db.newFluidIntFieldList(0, "surface point")
    vol = db.newFluidScalarFieldList(0.0, "volume")
    deltaCentroid = db.newFluidVectorFieldList(sph.Vector.zero, "delta centroid")
    A = db.newFluidScalarFieldList(0.0, "A")
    B = db.newFluidVectorFieldList(sph.Vector.zero, "B")
    C = db.newFluidTensorFieldList(sph.Tensor.zero, "B")
    gradA = db.newFluidVectorFieldList(sph.Vector.zero, "gradA")
    gradB = db.newFluidTensorFieldList(sph.Tensor.zero, "gradB")
    gradC = db.newFluidThirdRankTensorFieldList(sph.ThirdRankTensor.zero, "gradC")
    m0 = db.newFluidScalarFieldList(0.0, "m0")
    m1 = db.newFluidVectorFieldList(sph.Vector.zero, "m1")
    m2 = db.newFluidSymTensorFieldList(sph.SymTensor.zero, "m2")
    m3 = db.newFluidThirdRankTensorFieldList(sph.ThirdRankTensor.zero, "m3")
    m4 = db.newFluidFourthRankTensorFieldList(sph.FourthRankTensor.zero, "m4")
    gradm0 = db.newFluidVectorFieldList(sph.Vector.zero, "gradm0")
    gradm1 = db.newFluidTensorFieldList(sph.Tensor.zero, "gradm1")
    gradm2 = db.newFluidThirdRankTensorFieldList(sph.ThirdRankTensor.zero, "gradm2")
    gradm3 = db.newFluidFourthRankTensorFieldList(sph.FourthRankTensor.zero, "gradm3")
    gradm4 = db.newFluidFifthRankTensorFieldList(sph.FifthRankTensor.zero, "gradm4")

    if tessellationFileName is None:
        cells = sph.FacetedVolumeFieldList()
    else:
        cells = db.newFluidFacetedVolumeFieldList(sph.FacetedVolume(), "cells")

    # Kick start the volume using m/rho.
    for nodeListi, nodes in enumerate(db.fluidNodeLists()):
        for i in xrange(nodes.numInternalNodes):
            vol[nodeListi][i] = mass(nodeListi, i)/rhof(nodeListi, i)

    # We need the boundaries as a vector for calling iterateIdealH
    bound_vec = sph.vector_of_Boundary()
    for bc in boundaries:
        bound_vec.append(bc)

    # Iterate until we converge or max out.
    iter = 0
    avgdelta = 2.0*fracTol
    while (iter < 2) or (iter < maxIterations and avgdelta > fracTol) or mass.min() == 0.0:
        iter += 1

        # Remove any old ghost nodes info, and update the mass density
        for nodeListi, nodes in enumerate(db.fluidNodeLists()):
            nodes.numGhostNodes = 0
            nodes.neighbor().updateNodes()

            for i in xrange(nodes.numInternalNodes):
                rhof[nodeListi][i] = rhofunc(pos(nodeListi, i))

        # Create the new ghost nodes.
        for bc in boundaries:
            bc.setAllGhostNodes(db)
        for bc in boundaries:
            bc.finalizeGhostBoundary()
        for nodes in db.fluidNodeLists():
            nodes.neighbor().updateNodes()

        # Compute the new connectivity.
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()

        # Compute the new volumes and centroids (note this uses the old rho gradient, not quite right,
        # but expedient/efficient).
        sph.computeVoronoiVolume(pos, H, rhof, gradRhof, cm, W.kernelExtent, bounds, holes, surfacePoint, vol, deltaCentroid, cells)
        
        # Apply boundary conditions.
        for bc in boundaries:
            bc.applyFieldListGhostBoundary(vol)
            bc.applyFieldListGhostBoundary(rhof)
        for bc in boundaries:
            bc.finalizeGhostBoundary()

        # If the user provided a gradrho method, we can use it.  Otherwise we need to numerically evaluate
        # the density gradient.
        if gradrhofunc:
            for nodeListi, nodes in enumerate(db.fluidNodeLists()):
                for i in xrange(nodes.numInternalNodes):
                    gradRhof[nodeListi][i] = gradrhofunc(pos(nodeListi, i))

        else:
            # Use RK to numerically compute the new mass density gradient.
            sph.computeCRKSPHMoments(cm, W, vol, pos, H, correctionOrder, sph.NodeCoupling(),
                                     m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4)
            sph.computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, correctionOrder,
                                         A, B, C, gradA, gradB, gradC)
            gradRhof = sph.gradientCRKSPH(rhof, pos, vol, H, A, B, C, gradA, gradB, gradC, cm, correctionOrder, W)
        
        # Displace the points and update point masses.
        avgdelta = 0.0
        for nodeListi, nodes in enumerate(db.fluidNodeLists()):
            for i in xrange(nodes.numInternalNodes):
                delta = centroidFrac * deltaCentroid(nodeListi, i)
                if bounds:
                    while not bounds[nodeListi].contains(pos[nodeListi][i] + delta):
                        delta *= 0.9
                if holes:
                    for hole in holes[nodeListi]:
                        while hole.contains(pos[nodeListi][i] + delta):
                            delta *= 0.9
                if vol(nodeListi, i) > 0.0:
                    avgdelta += delta.magnitude()/vol(nodeListi, i)**(1.0/ndim)
                pos[nodeListi][i] += delta
                rhof[nodeListi][i] = rhofunc(pos(nodeListi, i))
                mass[nodeListi][i] = rhof(nodeListi,i)*vol(nodeListi,i)
        avgdelta = mpi.allreduce(avgdelta, mpi.SUM)/mpi.allreduce(db.numInternalNodes, mpi.SUM)
        print "centroidalRelaxNodes iteration %i, avg delta frac %g" % (iter, avgdelta)

        # Update the H tensors a bit.
        sph.iterateIdealH(db, bound_vec, W, sph.ASPHSmoothingScale(), 2)

    # Make a last pass updating the node data with the final info.
    for nodeListi, nodes in enumerate(db.fluidNodeLists()):
        nodes.numGhostNodes = 0
        nodes.neighbor().updateNodes()
        for i in xrange(nodes.numInternalNodes):
            rhof[nodeListi][i] = rhofunc(pos(nodeListi, i))
            mass[nodeListi][i] = rhof(nodeListi,i)*vol(nodeListi,i)

    # If requested, dump the final info to a diagnostic viz file.
    if tessellationFileName and SpheralVoronoiSiloDump:
        dumper = SpheralVoronoiSiloDump(baseFileName = tessellationFileName,
                                        listOfFieldLists = [vol, surfacePoint, mass, deltaCentroid],
                                        boundaries = boundaries,
                                        cells = cells)
        dumper.dump(0.0, iter)

    return vol, surfacePoint
