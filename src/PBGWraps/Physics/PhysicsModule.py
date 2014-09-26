from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

sys.path.append("../CXXTypes")
from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# Physics virtual bindings
#-------------------------------------------------------------------------------
def generatePhysicsVirtualBindings(x, ndim, pureVirtual):

    # Object names.
    vector = "Vector%id" % ndim
    tensor = "Tensor%id" % ndim
    symtensor = "SymTensor%id" % ndim
    scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
    database = "Spheral::DataBaseSpace::DataBase%id" % ndim
    state = "Spheral::State%id" % ndim
    derivatives = "Spheral::StateDerivatives%id" % ndim
    boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
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
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Physics/PhysicsTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("PhysicsSpace")

        # Expose types.
        self.Physics1d = addObject(space, "Physics1d", allow_subclassing=True)
        self.Physics2d = addObject(space, "Physics2d", allow_subclassing=True)
        self.Physics3d = addObject(space, "Physics3d", allow_subclassing=True)

        self.GenericHydro1d = addObject(space, "GenericHydro1d", allow_subclassing=True, parent=self.Physics1d)
        self.GenericHydro2d = addObject(space, "GenericHydro2d", allow_subclassing=True, parent=self.Physics2d)
        self.GenericHydro3d = addObject(space, "GenericHydro3d", allow_subclassing=True, parent=self.Physics3d)

        self.GenericBodyForce1d = addObject(space, "GenericBodyForce1d", allow_subclassing=True, parent=self.Physics1d)
        self.GenericBodyForce2d = addObject(space, "GenericBodyForce2d", allow_subclassing=True, parent=self.Physics2d)
        self.GenericBodyForce3d = addObject(space, "GenericBodyForce3d", allow_subclassing=True, parent=self.Physics3d)
        
        self.ArtificialConduction1d = addObject(space, "ArtificialConduction1d", allow_subclassing=True, parent=self.Physics1d)
        self.ArtificialConduction2d = addObject(space, "ArtificialConduction2d", allow_subclassing=True, parent=self.Physics2d)
        self.ArtificialConduction3d = addObject(space, "ArtificialConduction3d", allow_subclassing=True, parent=self.Physics3d)

        self.vector_of_Physics1d = addObject(mod, "vector_of_Physics1d", allow_subclassing=True)
        self.vector_of_Physics2d = addObject(mod, "vector_of_Physics2d", allow_subclassing=True)
        self.vector_of_Physics3d = addObject(mod, "vector_of_Physics3d", allow_subclassing=True)

        self.MassDensityType = space.add_enum("MassDensityType", ["SumDensity", "RigorousSumDensity", "HybridSumDensity", "IntegrateDensity", "VoronoiCellDensity", "SumVoronoiCellDensity"])
        self.HEvolutionType = space.add_enum("HEvolutionType", ["IdealH", "IntegrateH"])

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generatePhysicsBindings(self.Physics1d, 1)
        self.generatePhysicsBindings(self.Physics2d, 2)
        self.generatePhysicsBindings(self.Physics3d, 3)

        self.generateGenericHydroBindings(self.GenericHydro1d, 1)
        self.generateGenericHydroBindings(self.GenericHydro2d, 2)
        self.generateGenericHydroBindings(self.GenericHydro3d, 3)

        self.generateGenericBodyForceBindings(self.GenericBodyForce1d, 1)
        self.generateGenericBodyForceBindings(self.GenericBodyForce2d, 2)
        self.generateGenericBodyForceBindings(self.GenericBodyForce3d, 3)
        
        self.generateArtificialConductionBindings(self.ArtificialConduction1d, 1)
        self.generateArtificialConductionBindings(self.ArtificialConduction2d, 2)
        self.generateArtificialConductionBindings(self.ArtificialConduction3d, 3)

        generateStdVectorBindings(self.vector_of_Physics1d, "Spheral::PhysicsSpace::Physics1d*", "vector_of_Physics1d")
        generateStdVectorBindings(self.vector_of_Physics2d, "Spheral::PhysicsSpace::Physics2d*", "vector_of_Physics2d")
        generateStdVectorBindings(self.vector_of_Physics3d, "Spheral::PhysicsSpace::Physics3d*", "vector_of_Physics3d")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["PhysicsSpace"]

    #---------------------------------------------------------------------------
    # Physics bindings
    #---------------------------------------------------------------------------
    def generatePhysicsBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::Physics%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim

        # Constructor.
        x.add_constructor([])

        # Methods.
        x.add_method("appendBoundary", None, [refparam(boundary, "boundary")])
        x.add_method("prependBoundary", None, [refparam(boundary, "boundary")])
        x.add_method("clearBoundaries", None, [])
        x.add_method("haveBoundary", "bool", [constrefparam(boundary, "boundary")], is_const=True)
        const_ref_return_value(x, me, "%s::boundaryConditions" % me, vector_of_boundary, [], "boundaryConditions")
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
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
        x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               constrefparam(derivatives, "derivatives")], is_const=True, is_virtual=True)
        x.add_method("requireConnectivity", "bool", [], is_const=True, is_virtual=False)
        x.add_method("extraEnergy", "double", [], is_const=True, is_virtual=True)
        x.add_method("extraMomentum", vector, [], is_const=True, is_virtual=True)

        # The abstract interface.
        generatePhysicsVirtualBindings(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # GenericHydro
    #---------------------------------------------------------------------------
    def generateGenericHydroBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::GenericHydro%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim

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
        me = "Spheral::PhysicsSpace::GenericBodyForce%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim

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

    #---------------------------------------------------------------------------
    # Artificial Conduction
    #---------------------------------------------------------------------------
    def generateArtificialConductionBindings(self, x, ndim):
        
        # Object names.
        me = "Spheral::PhysicsSpace::ArtificialConduction%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        
        # Constructor.
        x.add_constructor([param("double", "arCondAlpha", default_value="0.5")])
        
        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        generatePhysicsVirtualBindings(x,ndim,pureVirtual=False)
        
        
        return

