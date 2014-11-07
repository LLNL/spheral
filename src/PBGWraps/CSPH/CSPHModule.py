import sys
from pybindgen import *

from PBGutils import *
from ref_return_value import *

sys.path.append("../Physics")
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class CSPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"CSPH/CSPHTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral.add_cpp_namespace("CSPHSpace")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        generichydro1d = findObject(PhysicsSpace, "GenericHydro1d")
        generichydro2d = findObject(PhysicsSpace, "GenericHydro2d")
        generichydro3d = findObject(PhysicsSpace, "GenericHydro3d")

        # Expose types.
        self.CSPHHydroBase1d = addObject(self.space, "CSPHHydroBase1d", allow_subclassing=True, parent=generichydro1d)
        self.CSPHHydroBase2d = addObject(self.space, "CSPHHydroBase2d", allow_subclassing=True, parent=generichydro2d)
        self.CSPHHydroBase3d = addObject(self.space, "CSPHHydroBase3d", allow_subclassing=True, parent=generichydro3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateCSPHHydroBaseBindings(self.CSPHHydroBase1d, 1)
        self.generateCSPHHydroBaseBindings(self.CSPHHydroBase2d, 2)
        self.generateCSPHHydroBaseBindings(self.CSPHHydroBase3d, 3)
        
        self.generateDimBindings(mod, 1)
        self.generateDimBindings(mod, 2)
        self.generateDimBindings(mod, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["CSPHSpace"]

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

        # Helper to compute the CSPH kernel.
        self.space.add_function("CSPHKernel", "double", [constrefparam(tablekernel, "W"),
                                                         constrefparam(vector, "rij"),
                                                         constrefparam(vector, "etai"),
                                                         param("double", "Hdeti"),
                                                         constrefparam(vector, "etaj"),
                                                         param("double", "Hdetj"),
                                                         param("double", "Ai"),
                                                         constrefparam(vector, "Bi")],
                                template_parameters = [dim],
                                custom_name = "CSPHKernel%id" % ndim,
                                docstring = "Evaluate the CSPH corrected kernel.")

        # Simultaneously evaluate the CSPH kernel and it's gradient.
        self.space.add_function("CSPHKernelAndGradient%id" % ndim, None, [constrefparam(tablekernel, "W"),
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
                                                                          Parameter.new("double*", "WCSPH", direction=Parameter.DIRECTION_OUT),
                                                                          Parameter.new("double*", "gradWSPH", direction=Parameter.DIRECTION_OUT),
                                                                          refparam(vector, "gradWCSPH")],
                                docstring = "Evaluate the CSPH corrected kernel and gradient simultaneously.")

        # CSPH sum density.
        self.space.add_function("computeCSPHSumMassDensity", None,
                                [constrefparam(connectivitymap, "connectivityMap"),
                                 constrefparam(tablekernel, "W"),
                                 constrefparam(vectorfieldlist, "position"),
                                 constrefparam(scalarfieldlist, "mass"),
                                 constrefparam(symtensorfieldlist, "H"),
                                 constrefparam(vector_of_boundary, "boundaries"),
                                 refparam(scalarfieldlist, "massDensity")],
                                template_parameters = [dim],
                                custom_name = "computeCSPHSumMassDensity%id" % ndim)

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

        # CSPH corrections.
        self.space.add_function("computeCSPHCorrections", None,
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
                                custom_name = "computeCSPHCorrections%id" % ndim)

        # CSPH Scalar interpolation.
        self.space.add_function("interpolateCSPH", scalarfieldlist,
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
                                custom_name = "interpolateCSPH")

        # CSPH Vector interpolation.
        self.space.add_function("interpolateCSPH", vectorfieldlist,
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
                                custom_name = "interpolateCSPH")

        # CSPH Tensor interpolation.
        self.space.add_function("interpolateCSPH", tensorfieldlist,
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
                                custom_name = "interpolateCSPH")

        # CSPH SymTensor interpolation.
        self.space.add_function("interpolateCSPH", symtensorfieldlist,
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
                                custom_name = "interpolateCSPH")

        # CSPH ThirdRankTensor interpolation.
        self.space.add_function("interpolateCSPH", thirdranktensorfieldlist,
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
                                custom_name = "interpolateCSPH")

        # CSPH Scalar gradient.
        self.space.add_function("gradientCSPH", vectorfieldlist,
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
                                custom_name = "gradientCSPH")

        # CSPH Vector gradient.
        self.space.add_function("gradientCSPH", tensorfieldlist,
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
                                custom_name = "gradientCSPH")

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

        # Compute the H scaled volume for each point.
        Spheral.add_function("computeHVolumes", None,
                             [param("double", "kernelExtent"),
                              constrefparam(symtensorfieldlist, "H"),
                              refparam(scalarfieldlist, "volume")],
                             docstring = "Compute the H scaled volume for each point.")

        return

    #---------------------------------------------------------------------------
    # Bindings (CSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateCSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::CSPHSpace::CSPHHydroBase%id" % ndim
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
                           param("double", "nTensile", default_value="4.0"),
                           param("int", "momentumConserving", default_value="true")])

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
        x.add_instance_attribute("momentumConserving", "bool", getter="momentumConserving", setter="momentumConserving")
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

    # #---------------------------------------------------------------------------
    # # Bindings (TotalHydro).
    # #---------------------------------------------------------------------------
    # def generateTotalHydroBindings(self, x, ndim):

    #     # Object names.
    #     me = "Spheral::PhysicsSpace::TotalHydro%id" % ndim
    #     dim = "Spheral::Dim<%i>" % ndim
    #     vector = "Vector%id" % ndim
    #     tensor = "Tensor%id" % ndim
    #     symtensor = "SymTensor%id" % ndim
    #     fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
    #     intfield = "Spheral::FieldSpace::IntField%id" % ndim
    #     scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
    #     vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
    #     vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
    #     tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
    #     thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
    #     vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
    #     vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
    #     vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
    #     symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
    #     intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
    #     scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
    #     vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
    #     vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
    #     tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
    #     symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
    #     thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
    #     vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
    #     vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
    #     vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
    #     nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
    #     state = "Spheral::State%id" % ndim
    #     derivatives = "Spheral::StateDerivatives%id" % ndim
    #     database = "Spheral::DataBaseSpace::DataBase%id" % ndim
    #     connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
    #     key = "pair_NodeList%id_string" % ndim
    #     vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
    #     tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
    #     artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
    #     fileio = "Spheral::FileIOSpace::FileIO"

    #     # Constructors.
    #     x.add_constructor([constrefparam(tablekernel, "W"),
    #                        constrefparam(tablekernel, "WPi"),
    #                        refparam(artificialviscosity, "Q"),
    #                        param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
    #                        param("double", "hmin", default_value="1.0e-100"),
    #                        param("double", "hmax", default_value="1.0e100"),
    #                        param("double", "hratiomin", default_value="0.1")])

    #     # Methods.
    #     x.add_method("evaluateDerivatives", None, [param("const double", "time"),
    #                                                param("const double", "dt"),
    #                                                constrefparam(database, "dataBase"),
    #                                                constrefparam(state, "state"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_const=True, is_virtual=True)
    #     x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
    #                                            refparam(state, "state"),
    #                                            constrefparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
    #     x.add_method("registerState", None, [refparam(database, "dataBase"),
    #                                          refparam(state, "state")],
    #                  is_virtual=True)
    #     x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
    #                                                refparam(derivatives, "derivatives")],
    #                  is_virtual=True)
    #     x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
    #     x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
    #     x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
    #     x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
    #                                      refparam("std::string", "pathName")],
    #                  is_const = True,
    #                  is_virtual = True)
    #     x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
    #                                         refparam("std::string", "pathName")],
    #                  is_virtual = True)
        
    #     # Attributes.
    #     x.add_instance_attribute("HEvolution", "HEvolutionType", getter="HEvolution", setter="HEvolution")
    #     x.add_instance_attribute("hmin", "double", getter="hmin", setter="hmin")
    #     x.add_instance_attribute("hmax", "double", getter="hmax", setter="hmax")
    #     x.add_instance_attribute("hratiomin", "double", getter="hratiomin", setter="hratiomin")
    #     x.add_instance_attribute("Hideal", symtensorfieldlist, getter="Hideal", is_const=True)
    #     x.add_instance_attribute("timeStepMask", intfieldlist, getter="timeStepMask", is_const=True)
    #     x.add_instance_attribute("pressure", scalarfieldlist, getter="pressure", is_const=True)
    #     x.add_instance_attribute("soundSpeed", scalarfieldlist, getter="soundSpeed", is_const=True)
    #     x.add_instance_attribute("weightedNeighborSum", scalarfieldlist, getter="weightedNeighborSum", is_const=True)
    #     x.add_instance_attribute("totalEnergy", scalarfieldlist, getter="totalEnergy", is_const=True)
    #     x.add_instance_attribute("volume", scalarfieldlist, getter="volume", is_const=True)
    #     x.add_instance_attribute("linearMomentum", vectorfieldlist, getter="linearMomentum", is_const=True)
    #     x.add_instance_attribute("DEDt", scalarfieldlist, getter="DEDt", is_const=True)
    #     x.add_instance_attribute("massSecondMoment", symtensorfieldlist, getter="massSecondMoment", is_const=True)

    #     return
