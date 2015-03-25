from SolidSpheral3d import *

#-------------------------------------------------------------------------------
# A package which keeps the H tensors eigen vector associated with hmin aligned
# with the local surface normal.
#-------------------------------------------------------------------------------
# class FacetedSurfaceASPHSmoothingScale(Physics):

#     def __init__(self, surface):
#         self.surface = surface
#         Physics.__init__(self)
#         return

#     def evaluateDerivatives(self, t, dt, db, state, derivs):
#         return

#     def dt(self, db, state, derivs, t):
#         return pair_double_string(1e100, "No vote")

#     def registerState(self, db, state):
#         return

#     def registerDerivatives(self, db, derivs):
#         return

#     def label(self):
#         return "FacetedSurfaceASPHSmoothingScale"

#     def finalize(self, t, dt, db, state, derivs):
#         pos = state.vectorFields(HydroFieldNames.position)
#         H = state.symTensorFields(HydroFieldNames.H)
#         for nodeListi in xrange(pos.numFields):
#             for i in xrange(pos[nodeListi].numInternalNodes):
#                 havg2 = (3.0/H(nodeListi, i).Trace())**2
#                 norm = Vector()
#                 for facet self.surface.facets():
#                     r2 = (facet.position - pos(nodeListi, i)).magnitude2()
#                     wi = exp(-r2/havg2)
#                     norm += facet.normal * wi
#                 norm = norm.unitVector()
#                 T = rotationMatrix(norm)
#                 Ti = T.Transpose()
#                 if 

#-------------------------------------------------------------------------------
# A special version of the ASPH smoothing scale algorithm designed to keep the 
# hmin eigen vector tensor aligned with the local faceted surface normal.
#-------------------------------------------------------------------------------
class FacetedSurfaceASPHSmoothingScale(ASPHSmoothingScale):
    def __init__(self,
                 surface,
                 nodes2facets):
        self.surface = surface
        self.nodes2facets = nodes2facets
        self.poshash2facet = {}
        ASPHSmoothingScale.__init__(self)
        return

    def idealSmoothingScale(self,
                            H,
                            pos,
                            zerothMoment,
                            secondMoment,
                            W,
                            hmin,
                            hmax,
                            hminratio,
                            nPerh,
                            connectivityMap,
                            nodeListi,
                            i):

        # Find the effective surface normal we want to use.
        fhat = Vector()
        verts = self.surface.vertices()
        facets = self.surface.facets()
        vertNorms = self.surface.vertexUnitNorms()
        hashi = hash(tuple(pos))
        assert hashi in self.poshash2facet
        facet = facets[self.poshash2facet[hashi]]
        # fhat = facet.normal/max(1.0e-10, (pos - facet.position).magnitude2())
        # for i in facet.ipoints:
        #     fhat += vertNorms[i]/max(1.0e-10, (pos - verts[i]).magnitude2())
        # fhat = fhat.unitVector()
        fhat = facet.normal
        T1 = rotationMatrix(fhat)
        
        # Now decompose the trial H tensor and align hmin with the surface normal.
        H0 = ASPHSmoothingScale.idealSmoothingScale(self,
                                                    H,
                                                    pos,
                                                    zerothMoment,
                                                    secondMoment,
                                                    W,
                                                    hmin,
                                                    hmax,
                                                    hminratio,
                                                    nPerh,
                                                    connectivityMap,
                                                    nodeListi,
                                                    i)
        Heigen = H0.eigenVectors()
        if ((Heigen.eigenValues.x > Heigen.eigenValues.y) and
            (Heigen.eigenValues.x > Heigen.eigenValues.z)):
            hvec = Heigen.eigenVectors.getColumn(0)
        elif ((Heigen.eigenValues.y > Heigen.eigenValues.x) and
              (Heigen.eigenValues.y > Heigen.eigenValues.z)):
            hvec = Heigen.eigenVectors.getColumn(1)
        else:
            hvec = Heigen.eigenVectors.getColumn(2)
        T0 = rotationMatrix(hvec)
        H0.rotationalTransform(T0)
        H0.rotationalTransform(T1.Transpose())
        return H0

#-------------------------------------------------------------------------------
# The specialized hydro object using our very particular method of H adaptation.
#-------------------------------------------------------------------------------
class FacetedSurfaceASPHHydro(SPHHydroBase):

    # Constructor.
    def __init__(self,
                 surface,
                 nodes2facets,
                 W,
                 WPi,
                 Q,
                 filter = 0.0,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = False,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.3,
                 nTensile = 4.0,
                 xmin = Vector(-1e100, -1e100, -1e100),
                 xmax = Vector( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = FacetedSurfaceASPHSmoothingScale(surface, nodes2facets)
        if xmin is None:
            xmin = Vector.one
        SPHHydroBase.__init__(self,
                              self._smoothingScaleMethod,
                              W,
                              WPi,
                              Q,
                              filter,
                              cfl,
                              useVelocityMagnitudeForDt,
                              compatibleEnergyEvolution,
                              gradhCorrection,
                              XSPH,
                              correctVelocityGradient,
                              sumMassDensityOverAllNodeLists,
                              densityUpdate,
                              HUpdate,
                              epsTensile,
                              nTensile,
                              xmin,
                              xmax)
        return

    # Override the pre-evaluateDerivative initialize step to update the 
    # node->facet mapping and build a lookup based on hashed position.  This
    # is necessary 'cause we don't pass the node index into the idealSmoothingScale 
    # algorithm above.
    def initialize(self, t, dt, dataBase, state, derivs):
        surface = self._smoothingScaleMethod.surface
        facets = surface.facets()
        nodes2facets = self._smoothingScaleMethod.nodes2facets
        facets2facets = surface.facetFacetConnectivity()

        # First update the node->facet info.  We assume here the node can't have moved
        # too much since the last update.
        poshash2facet = {}
        posfl = state.vectorFields(HydroFieldNames.position)
        assert len(posfl) == 1
        pos = posfl[0]
        n = pos.numInternalElements
        for i in xrange(n):
            rmin = 1e300
            fimin = None
            for fi in facets2facets[nodes2facets[i]]:
                r = facets[fi].distance(pos[i])
                if r < rmin:
                    rmin = r
                    fimin = fi
            nodes2facets[i] = fimin
            poshash2facet[hash(tuple(pos[i]))] = fimin

        self._smoothingScaleMethod.poshash2facet = poshash2facet
        SPHHydroBase.initialize(self, t, dt, dataBase, state, derivs)
        return

    # We need to restart some information in addition to the standard hydro variables.
    def label(self):
        return "FacetedSurfaceASPHHydro"

    def dumpState(self, file, path):
        file.write(self._smoothingScaleMethod.nodes2facets, path + "/nodes2facets")
        SPHHydroBase.dumpState(self, file, path)
        return

    def restoreState(self, file, path):
        file.read(self._smoothingScaleMethod.nodes2facets, path + "/nodes2facets")
        SPHHydroBase.restoreState(self, file, path)
        return

