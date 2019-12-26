import sys
from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class CRKSPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/CRKSPHTypes.hh"' % srcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")

        for dim in self.dims:
            exec('''
generichydro%(dim)id = findObject(self.space, "GenericHydro%(dim)id")
self.CRKSPHHydroBase%(dim)id = addObject(self.space, "CRKSPHHydroBase%(dim)id", allow_subclassing=True, parent=generichydro%(dim)id)
self.SolidCRKSPHHydroBase%(dim)id = addObject(self.space, "SolidCRKSPHHydroBase%(dim)id", allow_subclassing=True, parent=self.CRKSPHHydroBase%(dim)id)
self.CRKSPHVariant%(dim)id = addObject(self.space, "CRKSPHVariant%(dim)id", allow_subclassing=True, parent=self.CRKSPHHydroBase%(dim)id)
''' % {"dim" : dim})
        self.CRKOrder = self.space.add_enum("CRKOrder", [("ZerothOrder",    "Spheral::CRKOrder::ZerothOrder"),
                                                         ("LinearOrder",    "Spheral::CRKOrder::LinearOrder"),
                                                         ("QuadraticOrder", "Spheral::CRKOrder::QuadraticOrder")])
        self.CRKVolumeType = self.space.add_enum("CRKVolumeType", [("RKMassOverDensity", "Spheral::CRKVolumeType::RKMassOverDensity"),
                                                                   ("RKSumVolume",       "Spheral::CRKVolumeType::RKSumVolume"),
                                                                   ("RKVoronoiVolume",   "Spheral::CRKVolumeType::RKVoronoiVolume"),
                                                                   ("RKHullVolume",      "Spheral::CRKVolumeType::RKHullVolume"),
                                                                   ("HVolume",            "Spheral::CRKVolumeType::HVolume")])

        if 2 in self.dims:
            self.CRKSPHHydroBaseRZ = addObject(self.space, "CRKSPHHydroBaseRZ", allow_subclassing=True, parent=self.CRKSPHHydroBase2d)
            self.SolidCRKSPHHydroBaseRZ = addObject(self.space, "SolidCRKSPHHydroBaseRZ", allow_subclassing=True, parent=self.SolidCRKSPHHydroBase2d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generateCRKSPHHydroBaseBindings(self.CRKSPHHydroBase%(dim)id, %(dim)i)
self.generateSolidCRKSPHHydroBaseBindings(self.SolidCRKSPHHydroBase%(dim)id, %(dim)i)
self.generateCRKSPHVariantBindings(self.CRKSPHVariant%(dim)id, %(dim)i)
''' % {"dim" : dim})
            self.generateDimBindings(mod, dim)

        if 2 in self.dims:
            self.generateCRKSPHHydroBaseRZBindings()
            self.generateSolidCRKSPHHydroBaseRZBindings()

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, mod, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        thirdranktensor = "ThirdRankTensor%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        fourthranktensorfieldlist = "Spheral::FourthRankTensorFieldList%id" % ndim
        fifthranktensorfieldlist = "Spheral::FifthRankTensorFieldList%id" % ndim
        polyvolfieldlist = "Spheral::FacetedVolumeFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        vector_of_FacetedVolume = "vector_of_FacetedVolume%id" % ndim
        vector_of_vector_of_FacetedVolume = "vector_of_vector_of_FacetedVolume%id" % ndim
        polyvol = {1: "Box1d", 
                   2: "Polygon",
                   3: "Polyhedron"}[ndim]

        # Helper to compute the CRKSPH kernel.
        self.space.add_function("CRKSPHKernel", "double", [constrefparam(tablekernel, "W"),
                                                           param("Spheral::CRKOrder","correctionOrder"),
                                                           constrefparam(vector, "rij"),
                                                           constrefparam(vector, "etai"),
                                                           param("double", "Hdeti"),
                                                           constrefparam(vector, "etaj"),
                                                           param("double", "Hdetj"),
                                                           param("double", "Ai"),
                                                           constrefparam(vector, "Bi"),
                                                           constrefparam(tensor, "Ci"),
                                                           param("double", "correctionMin", default_value = "std::numeric_limits<double>::lowest()"),
                                                           param("double", "correctionMax", default_value = "std::numeric_limits<double>::max()")],
                                template_parameters = [dim],
                                custom_name = "CRKSPHKernel%id" % ndim,
                                docstring = "Evaluate the CRKSPH corrected kernel.")

        # Simultaneously evaluate the CRKSPH kernel and it's gradient.
        self.space.add_function("CRKSPHKernelAndGradient%id" % ndim, None, [Parameter.new("double*", "WCRKSPH", direction=Parameter.DIRECTION_OUT),
                                                                            Parameter.new("double*", "gradWSPH", direction=Parameter.DIRECTION_OUT),
                                                                            refparam(vector, "gradWCRKSPH"),
                                                                            constrefparam(tablekernel, "W"),
                                                                            param("Spheral::CRKOrder","correctionOrder"),
                                                                            constrefparam(vector, "rij"),
                                                                            constrefparam(vector, "etai"),
                                                                            constrefparam(symtensor, "Hi"),
                                                                            param("double", "Hdeti"),
                                                                            constrefparam(vector, "etaj"),
                                                                            constrefparam(symtensor, "Hj"),
                                                                            param("double", "Hdetj"),
                                                                            param("double", "Ai"),
                                                                            constrefparam(vector, "Bi"),
                                                                            constrefparam(tensor, "Ci"),
                                                                            constrefparam(vector, "gradAi"),
                                                                            constrefparam(tensor, "gradBi"),
                                                                            constrefparam(thirdranktensor, "gradCi"),
                                                                            param("double", "correctionMin", default_value = "std::numeric_limits<double>::lowest()"),
                                                                            param("double", "correctionMax", default_value = "std::numeric_limits<double>::max()")],
                                docstring = "Evaluate the CRKSPH corrected kernel and gradient simultaneously.")

        # CRKSPH sum density.
        self.space.add_function("computeCRKSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(scalarfieldlist, "vol"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(intfieldlist, "voidPoint"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHSumMassDensity%id" % ndim)

        self.space.add_function("computeVoronoiVolume", None,
                                [constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(scalarfieldlist, "rho"),
                                 constrefparam(vectorfieldlist, "gradRho"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(symtensorfieldlist, "damage"),
                                 constrefparam(vector_of_FacetedVolume, "facetedBoundaries"),
                                 constrefparam(vector_of_vector_of_FacetedVolume, "holes"),
                                 constrefparam(vector_of_boundary, "boundaries"),
                                 constrefparam(scalarfieldlist, "weights"),
                                 constrefparam(intfieldlist, "voidPoint"),
                                 refparam(intfieldlist, "surfacePoint"),
                                 refparam(scalarfieldlist, "vol"),
                                 refparam(vectorfieldlist, "deltaCentroid"),
                                 refparam(vectorvectorfieldlist, "etaVoidPoints"),
                                 refparam(polyvolfieldlist, "cells")],
                                docstring = "Compute the volume per point using a Voronoi tessellation.")
                                
        self.space.add_function("computeCRKSPHSumVolume", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "vol")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHSumVolume%id" % ndim)

        self.space.add_function("computeOccupancyVolume", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "vol")],
                                template_parameters = [dim],
                                custom_name = "computeOccupancyVolume%id" % ndim)

        self.space.add_function("computeSolidCRKSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(scalarfieldlist, "massDensity0"),
                                 constrefparam("Spheral::NodeCoupling", "nodeCoupling"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeSolidCRKSPHSumMassDensity%id" % ndim)

        # CRKSPH moments.
        self.space.add_function("computeCRKSPHMoments", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("Spheral::CRKOrder","correctionOrder"),
                                 constrefparam("Spheral::NodeCoupling", "nodeCoupling"),
                                 refparam(scalarfieldlist, "m0"),
                                 refparam(vectorfieldlist, "m1"),
                                 refparam(symtensorfieldlist, "m2"),
                                 refparam(thirdranktensorfieldlist, "m3"),
                                 refparam(fourthranktensorfieldlist, "m4"),
                                 refparam(vectorfieldlist, "gradm0"),
                                 refparam(tensorfieldlist, "gradm1"),
                                 refparam(thirdranktensorfieldlist, "gradm2"),
                                 refparam(fourthranktensorfieldlist, "gradm3"),
                                 refparam(fifthranktensorfieldlist, "gradm4")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHMoments%id" % ndim)
                                
        # Detect Surfaces
        self.space.add_function("detectSurface", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(scalarfieldlist, "m0"),
                                 constrefparam(vectorfieldlist, "m1"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam("double", "detectThreshold"),
                                 constrefparam("double", "detectRange"),
                                 constrefparam("double", "sweepAngle"),
                                 refparam(intfieldlist, "surfacePoint")],
                                template_parameters = [dim],
                                custom_name = "detectSurface%id" % ndim)

        # CRKSPH corrections.
        self.space.add_function("computeCRKSPHCorrections", None,
                                [constrefparam(scalarfieldlist, "m0"),
                                 constrefparam(vectorfieldlist, "m1"),
                                 constrefparam(symtensorfieldlist, "m2"),
                                 constrefparam(thirdranktensorfieldlist, "m3"),
                                 constrefparam(fourthranktensorfieldlist, "m4"),
                                 constrefparam(vectorfieldlist, "gradm0"),
                                 constrefparam(tensorfieldlist, "gradm1"),
                                 constrefparam(thirdranktensorfieldlist, "gradm2"),
                                 constrefparam(fourthranktensorfieldlist, "gradm3"),
                                 constrefparam(fifthranktensorfieldlist, "gradm4"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("Spheral::CRKOrder","correctionOrder"),
                                 refparam(scalarfieldlist, "A"),
                                 refparam(vectorfieldlist, "B"),
                                 refparam(tensorfieldlist, "C"),
                                 refparam(vectorfieldlist, "gradA"),
                                 refparam(tensorfieldlist, "gradB"),
                                 refparam(thirdranktensorfieldlist, "gradC")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHCorrections%id" % ndim)

        # CRKSPH interpolation.
        for (fl, element) in ((scalarfieldlist, "double"),
                              (vectorfieldlist, vector),
                              (tensorfieldlist, tensor),
                              (symtensorfieldlist, symtensor),
                              (thirdranktensorfieldlist, thirdranktensor)):
            self.space.add_function("interpolateCRKSPH", fl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(scalarfieldlist, "weight"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(scalarfieldlist, "A"),
                                     constrefparam(vectorfieldlist, "B"),
                                     constrefparam(tensorfieldlist, "C"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     param("Spheral::CRKOrder","correctionOrder"),
                                     constrefparam(tablekernel, "W")],
                                    template_parameters = [dim, element],
                                    custom_name = "interpolateCRKSPH",
                                    docstring = "interpolateCRKSPH: returns the CRK interpolations of the input FieldList (assumes full coupling between nodes).")
            self.space.add_function("interpolateCRKSPH", fl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(scalarfieldlist, "weight"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(scalarfieldlist, "A"),
                                     constrefparam(vectorfieldlist, "B"),
                                     constrefparam(tensorfieldlist, "C"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     param("Spheral::CRKOrder","correctionOrder"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam("Spheral::NodeCoupling", "nodeCoupling")],
                                    template_parameters = [dim, element],
                                    custom_name = "interpolateCRKSPH",
                                    docstring = "interpolateCRKSPH: returns the CRK interpolations of the input FieldList.")

        # CRKSPH gradient.
        for (fl, result, element) in ((scalarfieldlist, vectorfieldlist, "double"),
                                      (vectorfieldlist, tensorfieldlist, vector)):
            self.space.add_function("gradientCRKSPH", result,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(scalarfieldlist, "weight"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(scalarfieldlist, "A"),
                                     constrefparam(vectorfieldlist, "B"),
                                     constrefparam(tensorfieldlist, "C"),
                                     constrefparam(vectorfieldlist, "gradA"),
                                     constrefparam(tensorfieldlist, "gradB"),
                                     constrefparam(thirdranktensorfieldlist, "gradC"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     param("Spheral::CRKOrder","correctionOrder"),
                                     constrefparam(tablekernel, "W")],
                                    template_parameters = [dim, element],
                                    custom_name = "gradientCRKSPH",
                                    docstring = "Return the CRK gradient of the input FieldList (assumes full coupling between nodes).")
            self.space.add_function("gradientCRKSPH", result,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(scalarfieldlist, "weight"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(scalarfieldlist, "A"),
                                     constrefparam(vectorfieldlist, "B"),
                                     constrefparam(tensorfieldlist, "C"),
                                     constrefparam(vectorfieldlist, "gradA"),
                                     constrefparam(tensorfieldlist, "gradB"),
                                     constrefparam(thirdranktensorfieldlist, "gradC"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     param("Spheral::CRKOrder","correctionOrder"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam("Spheral::NodeCoupling", "nodeCoupling")],
                                    template_parameters = [dim, element],
                                    custom_name = "gradientCRKSPH",
                                    docstring = "Return the CRK gradient of the input FieldList.")

        # Center of mass with linear density gradient.
        self.space.add_function("centerOfMass", vector,
                                [constrefparam(polyvol, "polyvol"),
                                 constrefparam(vector, "gradRhoi")],
                                docstring = "Compute the center of mass for a %s assuming a linear density field." % polyvol)

        # Compute the hull volume for each point.
        Spheral = mod.add_cpp_namespace("Spheral")
        Spheral.add_function("computeHullVolumes", None,
                             [constrefparam(connectivitymap, "connectivityMap"),
                              param("double", "kernelExtent"),
                              constrefparam(vectorfieldlist, "position"),
                              constrefparam(symtensorfieldlist, "H"),
                              refparam(scalarfieldlist, "volume")],
                             docstring = "Compute the hull volume for each point in a FieldList of positions.")

        # Compute the H scaled volume for each point.
        Spheral.add_function("computeHVolumes", None,
                             [param("double", "nPerh"),
                              constrefparam(symtensorfieldlist, "H"),
                              refparam(scalarfieldlist, "volume")],
                             docstring = "Compute the H scaled volume for each point.")

        # Compute the hull volume for a neighbor set.
        Spheral.add_function("computeNeighborHull", polyvol,
                             [constrefparam("vector_of_vector_of_int", "fullConnectivity"),
                              param("double", "etaCutoff"),
                              constrefparam(vector, "ri"),
                              constrefparam(symtensor, "Hi"),
                              constrefparam(vectorfieldlist, "position")],
                             docstring = "Compute the hull volume for a given set of neighbors.")

        return

    #---------------------------------------------------------------------------
    # Bindings (CRKSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateCRKSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::CRKSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        voidboundary = "Spheral::CRKSPHVoidBoundary%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "evolveTotalEnergy"),
                           param("int", "XSPH"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("CRKOrder", "correctionOrder"),
                           param("CRKVolumeType", "volumeType"),
                           param("double", "epsTensile"),
                           param("double", "nTensile")])

        # Override the pure virtal overrides.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)

        # x.add_method("registerState", None, [refparam(database, "dataBase"),
        #                                      refparam(state, "state")],
        #              is_virtual=True)
        # x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
        #                                            refparam(derivatives, "derivatives")],
        #              is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        # x.add_method("evaluateDerivatives", None, [param("const double", "time"),
        #                                            param("const double", "dt"),
        #                                            constrefparam(database, "dataBase"),
        #                                            constrefparam(state, "state"),
        #                                            refparam(derivatives, "derivatives")],
        #              is_const=True, is_virtual=True)
        x.add_method("finalizeDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
        x.add_method("finalize", None, [param("const double", "time"),
                                        param("const double", "dt"),
                                        refparam(database, "dataBase"),
                                        refparam(state, "state"),
                                        refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        # x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        x.add_instance_attribute("densityUpdate", "MassDensityType", getter="densityUpdate", setter="densityUpdate")
        x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
        x.add_instance_attribute("correctionOrder", "CRKOrder", getter="correctionOrder", setter="correctionOrder")
        x.add_instance_attribute("volumeType", "CRKVolumeType", getter="volumeType", setter="volumeType")
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("evolveTotalEnergy", "bool", getter="evolveTotalEnergy", setter="evolveTotalEnergy")
        x.add_instance_attribute("XSPH", "bool", getter="XSPH", setter="XSPH")
        x.add_instance_attribute("filter", "double", getter="filter", setter="filter")
        x.add_instance_attribute("correctionMin", "double", getter="correctionMin", setter="correctionMin")
        x.add_instance_attribute("correctionMax", "double", getter="correctionMax", setter="correctionMax")

        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::entropy" % me, scalarfieldlist, [], "entropy")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
        const_ref_return_value(x, me, "%s::effectiveViscousPressure" % me, scalarfieldlist, [], "effectiveViscousPressure")
        const_ref_return_value(x, me, "%s::viscousWork" % me, scalarfieldlist, [], "viscousWork")
        const_ref_return_value(x, me, "%s::weightedNeighborSum" % me, scalarfieldlist, [], "weightedNeighborSum")
        const_ref_return_value(x, me, "%s::massSecondMoment" % me, symtensorfieldlist, [], "massSecondMoment")
        const_ref_return_value(x, me, "%s::volume" % me, scalarfieldlist, [], "volume")
        const_ref_return_value(x, me, "%s::massDensityGradient" % me, vectorfieldlist, [], "massDensityGradient")
        const_ref_return_value(x, me, "%s::XSPHDeltaV" % me, vectorfieldlist, [], "XSPHDeltaV")
        const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
        const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")
        const_ref_return_value(x, me, "%s::DmassDensityDt" % me, scalarfieldlist, [], "DmassDensityDt")
        const_ref_return_value(x, me, "%s::DspecificThermalEnergyDt" % me, scalarfieldlist, [], "DspecificThermalEnergyDt")
        const_ref_return_value(x, me, "%s::DHDt" % me, symtensorfieldlist, [], "DHDt")
        const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")
        const_ref_return_value(x, me, "%s::internalDvDx" % me, tensorfieldlist, [], "internalDvDx")
        const_ref_return_value(x, me, "%s::pairAccelerations" % me, vectorvectorfieldlist, [], "pairAccelerations")
        const_ref_return_value(x, me, "%s::deltaCentroid" % me, vectorfieldlist, [], "deltaCentroid")

        const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
        const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
        const_ref_return_value(x, me, "%s::C" % me, tensorfieldlist, [], "C")
        const_ref_return_value(x, me, "%s::gradA" % me, vectorfieldlist, [], "gradA")
        const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")
        const_ref_return_value(x, me, "%s::gradC" % me, thirdranktensorfieldlist, [], "gradC")
        const_ref_return_value(x, me, "%s::surfacePoint" % me, intfieldlist, [], "surfacePoint")
        const_ref_return_value(x, me, "%s::voidPoint" % me, intfieldlist, [], "voidPoint")
        const_ref_return_value(x, me, "%s::etaVoidPoints" % me, vectorvectorfieldlist, [], "etaVoidPoints")
        const_ref_return_value(x, me, "%s::m0" % me, scalarfieldlist, [], "m0")
        const_ref_return_value(x, me, "%s::m1" % me, vectorfieldlist, [], "m1")

        const_ref_return_value(x, me, "%s::voidBoundary" % me, voidboundary, [], "voidBoundary")

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidCRKSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateSolidCRKSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::SolidCRKSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "evolveTotalEnergy"),
                           param("int", "XSPH"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("CRKOrder", "correctionOrder"),
                           param("CRKVolumeType", "volumeType"),
                           param("double", "epsTensile"),
                           param("double", "nTensile"),
                           param("bool", "damageRelieveRubble", default_value="false")])

        # Attributes.
        const_ref_return_value(x, me, "%s::Hfield0" % me, symtensorfieldlist, [], "Hfield0")
        x.add_instance_attribute("damageRelieveRubble", "bool", getter="damageRelieveRubble", setter="damageRelieveRubble")

        return

    #---------------------------------------------------------------------------
    # Bindings (CRKSPHHydroBaseRZ).
    #---------------------------------------------------------------------------
    def generateCRKSPHHydroBaseRZBindings(self):

        # Object names.
        x = self.CRKSPHHydroBaseRZ
        ndim = 2
        me = "Spheral::CRKSPHHydroBaseRZ"
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "evolveTotalEnergy"),
                           param("int", "XSPH"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("CRKOrder", "correctionOrder"),
                           param("CRKVolumeType", "volumeType"),
                           param("double", "epsTensile"),
                           param("double", "nTensile")])

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidCRKSPHHydroBaseRZ).
    #---------------------------------------------------------------------------
    def generateSolidCRKSPHHydroBaseRZBindings(self):

        # Object names.
        x = self.SolidCRKSPHHydroBaseRZ
        ndim = 2
        me = "Spheral::SolidCRKSPHHydroBaseRZ"
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "evolveTotalEnergy"),
                           param("int", "XSPH"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("CRKOrder", "correctionOrder"),
                           param("CRKVolumeType", "volumeType"),
                           param("double", "epsTensile"),
                           param("double", "nTensile"),
                           param("bool", "damageRelieveRubble", default_value="false")])

        # Attributes.
        const_ref_return_value(x, me, "%s::deviatoricStressTT" % me, scalarfieldlist, [], "deviatoricStressTT")
        const_ref_return_value(x, me, "%s::DdeviatoricStressTTDt" % me, scalarfieldlist, [], "DdeviatoricStressTTDt")

        return

    #---------------------------------------------------------------------------
    # Bindings (CRKSPHVariant).
    #---------------------------------------------------------------------------
    def generateCRKSPHVariantBindings(self, x, ndim):

        # Object names.
        me = "Spheral::CRKSPHVariant%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIO"
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        voidboundary = "Spheral::CRKSPHVoidBoundary%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           param("double", "filter"),
                           param("double", "cfl"),
                           param("int", "useVelocityMagnitudeForDt"),
                           param("int", "compatibleEnergyEvolution"),
                           param("int", "evolveTotalEnergy"),
                           param("int", "XSPH"),
                           param("MassDensityType", "densityUpdate"),
                           param("HEvolutionType", "HUpdate"),
                           param("CRKOrder", "correctionOrder"),
                           param("CRKVolumeType", "volumeType"),
                           param("double", "epsTensile"),
                           param("double", "nTensile")])

        # Override the pure virtal overrides.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        return

