from pybindgen import *

from PBGutils import *
from ref_return_value import *

from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# Physics virtual bindings
#-------------------------------------------------------------------------------
def generatePhysicsVirtualBindings(x, ndim, pureVirtual):

    # Object names.
    vector = "Vector%id" % ndim
    tensor = "Tensor%id" % ndim
    symtensor = "SymTensor%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim
    database = "Spheral::DataBase%id" % ndim
    state = "Spheral::State%id" % ndim
    derivatives = "Spheral::StateDerivatives%id" % ndim
    boundary = "Spheral::Boundary%id" % ndim
    vector_of_boundary = "vector_of_Boundary%id" % ndim

    # Abstract Methods.
    x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                               param("const double", "dt"),
                                               constrefparam(database, "dataBase"),
                                               constrefparam(state, "state"),
                                               refparam(derivatives, "derivatives")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                              constrefparam(state, "state"),
                                              constrefparam(derivatives, "derivatives"),
                                              param("double", "time")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("registerState", None, [refparam(database, "dataBase"),
                                         refparam(state, "state")],
                 is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                               refparam(derivatives, "derivatives")],
                 is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("label", "std::string", [], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    return

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Physics:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/PhysicsTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        for dim in self.dims:
            exec('''
self.Physics%(dim)id = addObject(space, "Physics%(dim)id", allow_subclassing=True)
self.GenericHydro%(dim)id = addObject(space, "GenericHydro%(dim)id", allow_subclassing=True, parent=self.Physics%(dim)id)
self.GenericBodyForce%(dim)id = addObject(space, "GenericBodyForce%(dim)id", allow_subclassing=True, parent=self.Physics%(dim)id)
self.vector_of_Physics%(dim)id = addObject(mod, "vector_of_Physics%(dim)id", allow_subclassing=True)
''' % {"dim" : dim})

        self.MassDensityType = space.add_enum("MassDensityType", [("SumDensity", "Spheral::MassDensityType::SumDensity"),
                                                                  ("RigorousSumDensity", "Spheral::MassDensityType::RigorousSumDensity"),
                                                                  ("HybridSumDensity", "Spheral::MassDensityType::HybridSumDensity"),
                                                                  ("IntegrateDensity", "Spheral::MassDensityType::IntegrateDensity"),
                                                                  ("VoronoiCellDensity", "Spheral::MassDensityType::VoronoiCellDensity"),
                                                                  ("SumVoronoiCellDensity", "Spheral::MassDensityType::SumVoronoiCellDensity"),
                                                                  ("CorrectedSumDensity", "Spheral::MassDensityType::CorrectedSumDensity")])
        self.HEvolutionType = space.add_enum("HEvolutionType", [("IdealH", "Spheral::HEvolutionType::IdealH"),
                                                                ("IntegrateH", "Spheral::HEvolutionType::IntegrateH")])

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.generatePhysicsBindings(self.Physics%(dim)id, %(dim)i)
self.generateGenericHydroBindings(self.GenericHydro%(dim)id, %(dim)i)
self.generateGenericBodyForceBindings(self.GenericBodyForce%(dim)id, %(dim)i)
generateStdVectorBindings(self.vector_of_Physics%(dim)id, "Spheral::Physics%(dim)id*", "vector_of_Physics%(dim)id")
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Physics bindings
    #---------------------------------------------------------------------------
    def generatePhysicsBindings(self, x, ndim):

        # Object names.
        me = "Spheral::Physics%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim

        # Constructor.
        x.add_constructor([])

        # Methods.
        x.add_method("appendBoundary", None, [refparam(boundary, "boundary")])
        x.add_method("prependBoundary", None, [refparam(boundary, "boundary")])
        x.add_method("clearBoundaries", None, [])
        x.add_method("haveBoundary", "bool", [constrefparam(boundary, "boundary")], is_const=True)
        x.add_method("boundaryConditions", vector_of_boundary, [])
        #const_ref_return_value(x, me, "%s::boundaryConditions" % me, vector_of_boundary, [], "boundaryConditions")
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("preStepInitialize", None, [constrefparam(database, "dataBase"),
                                                 refparam(state, "state"),
                                                 refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("finalize", None, [param("const double", "time"),
                                        param("const double", "dt"),
                                        refparam(database, "dataBase"),
                                        refparam(state, "state"),
                                        refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("finalizeDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
        x.add_method("postStateUpdate", None, [param("double", "t"),
                                               param("double", "dt"),
                                               constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("requireConnectivity", "bool", [], is_const=True, is_virtual=False)
        x.add_method("requireGhostConnectivity", "bool", [], is_const=True, is_virtual=False)
        x.add_method("extraEnergy", "double", [], is_const=True, is_virtual=True)
        x.add_method("extraMomentum", vector, [], is_const=True, is_virtual=True)
        x.add_method("registerAdditionalVisualizationState", None, [refparam(database, "dataBase"),
                                                                    refparam(state, "state")],
                     is_virtual=True)

        # The abstract interface.
        generatePhysicsVirtualBindings(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # GenericHydro
    #---------------------------------------------------------------------------
    def generateGenericHydroBindings(self, x, ndim):

        # Object names.
        me = "Spheral::GenericHydro%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim

        # Constructor.
        x.add_constructor([constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           refparam(artificialviscosity, "Q"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false")])

        # Methods.
        x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                                  constrefparam(state, "state"),
                                                  constrefparam(derivatives, "derivatives"),
                                                  param("double", "time")],
                     is_const=True, is_virtual=True)
        ref_return_value(x, me, "%s::artificialViscosity" % me, artificialviscosity, [], "artificialViscosity")
        const_ref_return_value(x, me, "%s::kernel" % me, tablekernel, [], "kernel")
        const_ref_return_value(x, me, "%s::PiKernel" % me, tablekernel, [], "PiKernel")

        x.add_method("updateMasterNeighborStats", None, [param("int", "numMaster")], is_const=True, visibility="protected")
        x.add_method("updateCoarseNeighborStats", None, [param("int", "numCoarse")], is_const=True, visibility="protected")
        x.add_method("updateRefineNeighborStats", None, [param("int", "numRefine")], is_const=True, visibility="protected")
        x.add_method("updateActualNeighborStats", None, [param("int", "numActual")], is_const=True, visibility="protected")

        # Attributes.
        x.add_instance_attribute("cfl", "double", getter="cfl", setter="cfl")
        x.add_instance_attribute("useVelocityMagnitudeForDt", "bool", getter="useVelocityMagnitudeForDt", setter="useVelocityMagnitudeForDt")

        x.add_instance_attribute("minMasterNeighbor", "int", getter="minMasterNeighbor", is_const=True)
        x.add_instance_attribute("maxMasterNeighbor", "int", getter="maxMasterNeighbor", is_const=True)
        x.add_instance_attribute("averageMasterNeighbor", "double", getter="averageMasterNeighbor", is_const=True)

        x.add_instance_attribute("minCoarseNeighbor", "int", getter="minCoarseNeighbor", is_const=True)
        x.add_instance_attribute("maxCoarseNeighbor", "int", getter="maxCoarseNeighbor", is_const=True)
        x.add_instance_attribute("averageCoarseNeighbor", "double", getter="averageCoarseNeighbor", is_const=True)

        x.add_instance_attribute("minRefineNeighbor", "int", getter="minRefineNeighbor", is_const=True)
        x.add_instance_attribute("maxRefineNeighbor", "int", getter="maxRefineNeighbor", is_const=True)
        x.add_instance_attribute("averageRefineNeighbor", "double", getter="averageRefineNeighbor", is_const=True)

        x.add_instance_attribute("minActualNeighbor", "int", getter="minActualNeighbor", is_const=True)
        x.add_instance_attribute("maxActualNeighbor", "int", getter="maxActualNeighbor", is_const=True)
        x.add_instance_attribute("averageActualNeighbor", "double", getter="averageActualNeighbor", is_const=True)

        return

    #---------------------------------------------------------------------------
    # GenericBodyForce
    #---------------------------------------------------------------------------
    def generateGenericBodyForceBindings(self, x, ndim):

        # Object names.
        me = "Spheral::GenericBodyForce%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim

        # Constructor.
        x.add_constructor([])

        # Methods.
        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
        const_ref_return_value(x, me, "%s::DxDt" % me, vectorfieldlist, [], "DxDt")
        const_ref_return_value(x, me, "%s::DvDt" % me, vectorfieldlist, [], "DvDt")

        return
