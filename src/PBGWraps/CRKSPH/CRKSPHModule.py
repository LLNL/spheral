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
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral.add_cpp_namespace("CRKSPHSpace")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")

        for dim in self.dims:
            exec('''
generichydro%(dim)id = findObject(PhysicsSpace, "GenericHydro%(dim)id")
self.CRKSPHHydroBase%(dim)id = addObject(self.space, "CRKSPHHydroBase%(dim)id", allow_subclassing=True, parent=generichydro%(dim)id)
self.SolidCRKSPHHydroBase%(dim)id = addObject(self.space, "SolidCRKSPHHydroBase%(dim)id", allow_subclassing=True, parent=self.CRKSPHHydroBase%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generateCRKSPHHydroBaseBindings(self.CRKSPHHydroBase%(dim)id, %(dim)i)
self.generateSolidCRKSPHHydroBaseBindings(self.SolidCRKSPHHydroBase%(dim)id, %(dim)i)
''' % {"dim" : dim})
            self.generateDimBindings(mod, dim)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["CRKSPHSpace"]

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, mod, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        thirdranktensor = "ThirdRankTensor%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        polyvolfieldlist = "Spheral::FieldSpace::FacetedVolumeFieldList%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        polyvol = {1: "Box1d", 
                   2: "Polygon",
                   3: "Polyhedron"}[ndim]

        # Helper to compute the CRKSPH kernel.
        self.space.add_function("CRKSPHKernel", "double", [constrefparam(tablekernel, "W"),
                                                           constrefparam(vector, "rij"),
                                                           constrefparam(vector, "etai"),
                                                           param("double", "Hdeti"),
                                                           constrefparam(vector, "etaj"),
                                                           param("double", "Hdetj"),
                                                           param("double", "Ai"),
                                                           constrefparam(vector, "Bi")],
                                template_parameters = [dim],
                                custom_name = "CRKSPHKernel%id" % ndim,
                                docstring = "Evaluate the CRKSPH corrected kernel.")

        # Simultaneously evaluate the CRKSPH kernel and it's gradient.
        self.space.add_function("CRKSPHKernelAndGradient%id" % ndim, None, [constrefparam(tablekernel, "W"),
                                                                            constrefparam(vector, "rij"),
                                                                            constrefparam(vector, "etai"),
                                                                            constrefparam(symtensor, "Hi"),
                                                                            param("double", "Hdeti"),
                                                                            constrefparam(vector, "etaj"),
                                                                            constrefparam(symtensor, "Hj"),
                                                                            param("double", "Hdetj"),
                                                                            param("double", "Ai"),
                                                                            constrefparam(vector, "Bi"),
                                                                            constrefparam(vector, "gradAi"),
                                                                            constrefparam(tensor, "gradBi"),
                                                                            Parameter.new("double*", "WCRKSPH", direction=Parameter.DIRECTION_OUT),
                                                                            Parameter.new("double*", "gradWSPH", direction=Parameter.DIRECTION_OUT),
                                                                            refparam(vector, "gradWCRKSPH")],
                                docstring = "Evaluate the CRKSPH corrected kernel and gradient simultaneously.")

        # CRKSPH sum density.
        self.space.add_function("computeCRKSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHSumMassDensity%id" % ndim)

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

        # Hull sum density.
        self.space.add_function("computeHullSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam("Spheral::NodeCoupling", "nodeCoupling"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeHullSumMassDensity%id" % ndim)

        # CRKSPH corrections.
        self.space.add_function("computeCRKSPHCorrections", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(tablekernel, "Wf"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "A"),
                                 refparam(vectorfieldlist, "B"),
                                 refparam(vectorfieldlist, "gradA"),
                                 refparam(tensorfieldlist, "gradB"),
                                 refparam(scalarfieldlist, "Af"),
                                 refparam(vectorfieldlist, "Bf"),
                                 refparam(vectorfieldlist, "gradAf"),
                                 refparam(tensorfieldlist, "gradBf")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHCorrections%id" % ndim)
        self.space.add_function("computeCRKSPHCorrections", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam("Spheral::NodeCoupling", "nodeCoupling"),
                                 refparam(scalarfieldlist, "A"),
                                 refparam(vectorfieldlist, "B"),
                                 refparam(vectorfieldlist, "gradA"),
                                 refparam(tensorfieldlist, "gradB"),
                                 refparam(scalarfieldlist, "Ac"),
                                 refparam(vectorfieldlist, "Bc"),
                                 refparam(vectorfieldlist, "gradAc"),
                                 refparam(tensorfieldlist, "gradBc")],
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
                                     constrefparam(connectivitymap, "connectivityMap"),
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
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam("Spheral::NodeCoupling", "nodeCoupling")],
                                    template_parameters = [dim, element],
                                    custom_name = "interpolateCRKSPH",
                                    docstring = "interpolateCRKSPH: returns the CRK interpolations of the input FieldList.")

        # CRDSPH gradient.
        for (fl, result, element) in ((scalarfieldlist, vectorfieldlist, "double"),
                                      (vectorfieldlist, tensorfieldlist, vector)):
            self.space.add_function("gradientCRKSPH", result,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(scalarfieldlist, "weight"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(scalarfieldlist, "A"),
                                     constrefparam(vectorfieldlist, "B"),
                                     constrefparam(vectorfieldlist, "gradA"),
                                     constrefparam(tensorfieldlist, "gradB"),
                                     constrefparam(connectivitymap, "connectivityMap"),
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
                                     constrefparam(vectorfieldlist, "gradA"),
                                     constrefparam(tensorfieldlist, "gradB"),
                                     constrefparam(connectivitymap, "connectivityMap"),
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
                              refparam(polyvolfieldlist, "polyvol"),
                              refparam(scalarfieldlist, "volume")],
                             docstring = "Compute the hull volume for each point in a FieldList of positions.")

        # Compute the centroids.
        Spheral.add_function("computeVoronoiCentroids", vectorfieldlist,
                             [constrefparam(vectorfieldlist, "position")],
                             docstring = "Compute the Voronoi based centroids for each point in a FieldList of positions.")

        # Compute the H scaled volume for each point.
        Spheral.add_function("computeHVolumes", None,
                             [param("double", "kernelExtent"),
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
        me = "Spheral::CRKSPHSpace::CRKSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "Wfilter"),
                           param("double", "filter", default_value="0.1"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "XSPH", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0")])

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
        x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               constrefparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
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
        x.add_instance_attribute("compatibleEnergyEvolution", "bool", getter="compatibleEnergyEvolution", setter="compatibleEnergyEvolution")
        x.add_instance_attribute("XSPH", "bool", getter="XSPH", setter="XSPH")
        x.add_instance_attribute("filter", "double", getter="filter", setter="filter")

        const_ref_return_value(x, me, "%s::filterKernel" % me, tablekernel, [], "filterKernel")
        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::volume" % me, scalarfieldlist, [], "volume")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
        const_ref_return_value(x, me, "%s::effectiveViscousPressure" % me, scalarfieldlist, [], "effectiveViscousPressure")
        const_ref_return_value(x, me, "%s::viscousWork" % me, scalarfieldlist, [], "viscousWork")
        const_ref_return_value(x, me, "%s::weightedNeighborSum" % me, scalarfieldlist, [], "weightedNeighborSum")
        const_ref_return_value(x, me, "%s::massSecondMoment" % me, symtensorfieldlist, [], "massSecondMoment")
        const_ref_return_value(x, me, "%s::XSPHDeltaV" % me, vectorfieldlist, [], "XSPHDeltaV")
        const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
        const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")
        const_ref_return_value(x, me, "%s::DmassDensityDt" % me, scalarfieldlist, [], "DmassDensityDt")
        const_ref_return_value(x, me, "%s::DspecificThermalEnergyDt" % me, scalarfieldlist, [], "DspecificThermalEnergyDt")
        const_ref_return_value(x, me, "%s::DHDt" % me, symtensorfieldlist, [], "DHDt")
        const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")
        const_ref_return_value(x, me, "%s::internalDvDx" % me, tensorfieldlist, [], "internalDvDx")
        const_ref_return_value(x, me, "%s::DmassDensityDx" % me, vectorfieldlist, [], "DmassDensityDx")
        const_ref_return_value(x, me, "%s::pairAccelerations" % me, vectorvectorfieldlist, [], "pairAccelerations")

        const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
        const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
        const_ref_return_value(x, me, "%s::gradA" % me, vectorfieldlist, [], "gradA")
        const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")
        const_ref_return_value(x, me, "%s::surfNorm" % me, vectorfieldlist, [], "surfNorm")

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidCRKSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateSolidCRKSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::CRKSPHSpace::SolidCRKSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           refparam(artificialviscosity, "Q"),
                           constrefparam(tablekernel, "Wfilter"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "XSPH", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0")])

        # Attributes.
        const_ref_return_value(x, me, "%s::Adamage" % me, scalarfieldlist, [], "Adamage")
        const_ref_return_value(x, me, "%s::Bdamage" % me, vectorfieldlist, [], "Bdamage")
        const_ref_return_value(x, me, "%s::gradAdamage" % me, vectorfieldlist, [], "gradAdamage")
        const_ref_return_value(x, me, "%s::gradBdamage" % me, tensorfieldlist, [], "gradBdamage")

        return
