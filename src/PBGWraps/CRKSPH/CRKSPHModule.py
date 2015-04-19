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
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/CRKSPHTypes.hh"' % srcdir)
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral.add_cpp_namespace("CRKSPHSpace")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        generichydro1d = findObject(PhysicsSpace, "GenericHydro1d")
        generichydro2d = findObject(PhysicsSpace, "GenericHydro2d")
        generichydro3d = findObject(PhysicsSpace, "GenericHydro3d")

        # Expose types.
        self.CRKSPHHydroBase1d = addObject(self.space, "CRKSPHHydroBase1d", allow_subclassing=True, parent=generichydro1d)
        self.CRKSPHHydroBase2d = addObject(self.space, "CRKSPHHydroBase2d", allow_subclassing=True, parent=generichydro2d)
        self.CRKSPHHydroBase3d = addObject(self.space, "CRKSPHHydroBase3d", allow_subclassing=True, parent=generichydro3d)

        # Expose types.
        self.SolidCRKSPHHydroBase1d = addObject(self.space, "SolidCRKSPHHydroBase1d", allow_subclassing=True, parent=self.CRKSPHHydroBase1d)
        self.SolidCRKSPHHydroBase2d = addObject(self.space, "SolidCRKSPHHydroBase2d", allow_subclassing=True, parent=self.CRKSPHHydroBase2d)
        self.SolidCRKSPHHydroBase3d = addObject(self.space, "SolidCRKSPHHydroBase3d", allow_subclassing=True, parent=self.CRKSPHHydroBase3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateCRKSPHHydroBaseBindings(self.CRKSPHHydroBase1d, 1)
        self.generateCRKSPHHydroBaseBindings(self.CRKSPHHydroBase2d, 2)
        self.generateCRKSPHHydroBaseBindings(self.CRKSPHHydroBase3d, 3)
        
        self.generateDimBindings(mod, 1)
        self.generateDimBindings(mod, 2)
        self.generateDimBindings(mod, 3)

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
                                 constrefparam(vector_of_boundary, "boundaries"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHSumMassDensity%id" % ndim)

        # Hull sum density.
        self.space.add_function("computeHullSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeHullSumMassDensity%id" % ndim)

        # CRKSPH corrections.
        self.space.add_function("computeCRKSPHCorrections", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 refparam(scalarfieldlist, "m0"),
                                 refparam(vectorfieldlist, "m1"),
                                 refparam(symtensorfieldlist, "m2"),
                                 refparam(scalarfieldlist, "A0"),
                                 refparam(scalarfieldlist, "A"),
                                 refparam(vectorfieldlist, "B"),
                                 refparam(vectorfieldlist, "C"),
                                 refparam(tensorfieldlist, "D"),
                                 refparam(vectorfieldlist, "gradA0"),
                                 refparam(vectorfieldlist, "gradA"),
                                 refparam(tensorfieldlist, "gradB")],
                                template_parameters = [dim],
                                custom_name = "computeCRKSPHCorrections%id" % ndim)

        # CRKSPH Scalar interpolation.
        self.space.add_function("interpolateCRKSPH", scalarfieldlist,
                                [constrefparam(scalarfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, "double"],
                                custom_name = "interpolateCRKSPH")

        # CRKSPH Vector interpolation.
        self.space.add_function("interpolateCRKSPH", vectorfieldlist,
                                [constrefparam(vectorfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, vector],
                                custom_name = "interpolateCRKSPH")

        # CRKSPH Tensor interpolation.
        self.space.add_function("interpolateCRKSPH", tensorfieldlist,
                                [constrefparam(tensorfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, tensor],
                                custom_name = "interpolateCRKSPH")

        # CRKSPH SymTensor interpolation.
        self.space.add_function("interpolateCRKSPH", symtensorfieldlist,
                                [constrefparam(symtensorfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, symtensor],
                                custom_name = "interpolateCRKSPH")

        # CRKSPH ThirdRankTensor interpolation.
        self.space.add_function("interpolateCRKSPH", thirdranktensorfieldlist,
                                [constrefparam(thirdranktensorfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 param("bool", "coupleNodeLists"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, thirdranktensor],
                                custom_name = "interpolateCRKSPH")

        # CRKSPH Scalar gradient.
        self.space.add_function("gradientCRKSPH", vectorfieldlist,
                                [constrefparam(scalarfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(vectorfieldlist, "C"),
                                 constrefparam(tensorfieldlist, "D"),
                                 constrefparam(vectorfieldlist, "gradA"),
                                 constrefparam(tensorfieldlist, "gradB"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, "double"],
                                custom_name = "gradientCRKSPH")

        # CRKSPH Vector gradient.
        self.space.add_function("gradientCRKSPH", tensorfieldlist,
                                [constrefparam(vectorfieldlist, "fieldList"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "weight"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(scalarfieldlist, "A"),
                                 constrefparam(vectorfieldlist, "B"),
                                 constrefparam(vectorfieldlist, "C"),
                                 constrefparam(tensorfieldlist, "D"),
                                 constrefparam(vectorfieldlist, "gradA"),
                                 constrefparam(tensorfieldlist, "gradB"),
                                 constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W")],
                                template_parameters = [dim, vector],
                                custom_name = "gradientCRKSPH")

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
                           param("double", "filter", default_value="0.1"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "XSPH", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0")])

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
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
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
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

        const_ref_return_value(x, me, "%s::smoothingScaleMethod" % me, smoothingscalebase, [], "smoothingScaleMethod")
        const_ref_return_value(x, me, "%s::timeStepMask" % me, intfieldlist, [], "timeStepMask")
        const_ref_return_value(x, me, "%s::pressure" % me, scalarfieldlist, [], "pressure")
        const_ref_return_value(x, me, "%s::soundSpeed" % me, scalarfieldlist, [], "soundSpeed")
        const_ref_return_value(x, me, "%s::volume" % me, scalarfieldlist, [], "volume")
        const_ref_return_value(x, me, "%s::specificThermalEnergy0" % me, scalarfieldlist, [], "specificThermalEnergy0")
        const_ref_return_value(x, me, "%s::Hideal" % me, symtensorfieldlist, [], "Hideal")
        const_ref_return_value(x, me, "%s::maxViscousPressure" % me, scalarfieldlist, [], "maxViscousPressure")
        #const_ref_return_value(x, me, "%s::massDensitySum" % me, scalarfieldlist, [], "massDensitySum")
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

        const_ref_return_value(x, me, "%s::m0" % me, scalarfieldlist, [], "m0")
        const_ref_return_value(x, me, "%s::m1" % me, vectorfieldlist, [], "m1")
        const_ref_return_value(x, me, "%s::m2" % me, symtensorfieldlist, [], "m2")
        const_ref_return_value(x, me, "%s::A0" % me, scalarfieldlist, [], "A0")
        const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
        const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
        const_ref_return_value(x, me, "%s::C" % me, vectorfieldlist, [], "C")
        const_ref_return_value(x, me, "%s::D" % me, tensorfieldlist, [], "D")
        const_ref_return_value(x, me, "%s::gradA" % me, vectorfieldlist, [], "gradA")
        const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")

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
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.25"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "XSPH", default_value="true"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.0"),
                           param("double", "nTensile", default_value="4.0")])

        return
