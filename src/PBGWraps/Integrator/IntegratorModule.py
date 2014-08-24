from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Integrator:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Integrator/IntegratorTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("IntegratorSpace")

        # Expose types.
        for dim in xrange(3):
            exec('''
self.Integrator%(dim)id = addObject(space, "Integrator%(dim)id", allow_subclassing=True)
self.PredictorCorrector%(dim)id = addObject(space, "PredictorCorrectorIntegrator%(dim)id", parent=self.Integrator%(dim)id)
self.SynchronousRK1Integrator%(dim)id = addObject(space, "SynchronousRK1Integrator%(dim)id", parent=self.Integrator%(dim)id)
self.SynchronousRK2Integrator%(dim)id = addObject(space, "SynchronousRK2Integrator%(dim)id", parent=self.Integrator%(dim)id)
self.SynchronousRK4Integrator%(dim)id = addObject(space, "SynchronousRK4Integrator%(dim)id", parent=self.Integrator%(dim)id)
self.CheapSynchronousRK2Integrator%(dim)id = addObject(space, "CheapSynchronousRK2Integrator%(dim)id", parent=self.Integrator%(dim)id)
self.VerletIntegrator%(dim)id = addObject(space, "VerletIntegrator%(dim)id", parent=self.Integrator%(dim)id)
''' % {"dim" : dim + 1})

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in xrange(3):
            exec('''
self.generateIntegratorBindings(self.Integrator%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.PredictorCorrector%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.SynchronousRK1Integrator%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.SynchronousRK2Integrator%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.SynchronousRK4Integrator%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.CheapSynchronousRK2Integrator%(dim)id, %(dim)i)
self.generateIntegratorDescendentBindings(self.VerletIntegrator%(dim)id, %(dim)i)
''' % {"dim" : dim + 1})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["IntegratorSpace"]

    #---------------------------------------------------------------------------
    # Generate bindings (Integrator).
    #---------------------------------------------------------------------------
    def generateIntegratorBindings(self, x, ndim):

        # Object names.
        me = "Spheral::IntegratorSpace::Integrator%id" % ndim
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
        physics = "Spheral::PhysicsSpace::Physics%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        vectorphysics = "vector_of_Physics%id" % ndim
        vectorboundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(database, "dataBase")])
        x.add_constructor([refparam(database, "dataBase"),
                           constrefparam(vectorphysics, "physicsPackages")])

        # Methods.
        x.add_method("step", None, [param("double", "maxTime")], is_virtual=True, is_pure_virtual=True)
        x.add_method("selectDt", "double", [param("double", "dtMin"),
                                            param("double", "dtMax"),
                                            constrefparam(state, "state"),
                                            constrefparam(derivatives, "derivs")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("initialize", None, [refparam(state, "state"),
                                          refparam(derivatives, "derivs")],
                     is_virtual = True)
        x.add_method("preStepInitialize", None, [param("double", "t"),
                                                 param("double", "dt"),
                                                 refparam(state, "state"),
                                                 refparam(derivatives, "derivs")],
                     is_virtual = True)
        x.add_method("finalize", None, [param("double", "t"),
                                        param("double", "dt"),
                                        refparam(state, "state"),
                                        refparam(derivatives, "derivs")],
                     is_virtual = True)
        x.add_method("evaluateDerivatives", None, [param("double", "t"),
                                                   param("double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivs")],
                     is_const = True)
        x.add_method("finalizeDerivatives", None, [param("double", "t"),
                                                   param("double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivs")],
                     is_const = True)
        x.add_method("postStateUpdate", None, [constrefparam(database, "dataBase"),
                                               refparam(state, "state"),
                                               constrefparam(derivatives, "derivs")],
                     is_const = True)
        x.add_method("appendPhysicsPackage", None, [constrefparam(physics, "package")])
        x.add_method("havePhysicsPackage", "bool", [constrefparam(physics, "package")], is_const=True)
        x.add_method("uniqueBoundaryConditions", vectorboundary, [], is_const=True)
        x.add_method("setGhostNodes", None, [])
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivs")])
        x.add_method("finalizeGhostBoundaries", None, [])
        x.add_method("setViolationNodes", None, [])
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivs")])
        x.add_method("copyGhostState", None, [constrefparam(state, "state0"), refparam(state, "state1")], is_const=True)
        x.add_method("advance", None, [param("double", "goalTime")], is_virtual=True)
        const_ref_return_value(x, me, "%s::dataBase" % me, database, [], "dataBase")
        const_ref_return_value(x, me, "%s::physicsPackages" % me, vectorphysics, [], "physicsPackages")
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("currentTime", "double", getter="currentTime", setter="currentTime")
        x.add_instance_attribute("currentCycle", "int", getter="currentCycle", setter="currentCycle")
        x.add_instance_attribute("dtMin", "double", getter="dtMin", setter="dtMin")
        x.add_instance_attribute("dtMax", "double", getter="dtMax", setter="dtMax")
        x.add_instance_attribute("lastDt", "double", getter="lastDt", setter="lastDt")
        x.add_instance_attribute("dtGrowth", "double", getter="dtGrowth", setter="dtGrowth")
        x.add_instance_attribute("rigorousBoundaries", "bool", getter="rigorousBoundaries", setter="rigorousBoundaries")
        x.add_instance_attribute("updateBoundaryFrequency", "int", getter="updateBoundaryFrequency", setter="updateBoundaryFrequency")
        x.add_instance_attribute("verbose", "bool", getter="verbose", setter="verbose")
        x.add_instance_attribute("dtThreshold", "double", getter="dtThreshold", setter="dtThreshold")
        x.add_instance_attribute("domainDecompositionIndependent", "bool", getter="domainDecompositionIndependent", setter="domainDecompositionIndependent")
        x.add_instance_attribute("cullGhostNodes", "bool", getter="cullGhostNodes", setter="cullGhostNodes")

        return

    #---------------------------------------------------------------------------
    # Generate bindings (all descendents of Integrator).
    #---------------------------------------------------------------------------
    def generateIntegratorDescendentBindings(self, x, ndim):

        # Object names.
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        vectorphysics = "vector_of_Physics%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(database, "dataBase")])
        x.add_constructor([refparam(database, "dataBase"),
                           constrefparam(vectorphysics, "physicsPackages")])

        # You've got to provide this!
        x.add_method("step", None, [param("double", "maxTime")], is_virtual=True)
