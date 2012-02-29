from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Gravity:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Gravity/GravityTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        space = Spheral.add_cpp_namespace("GravitySpace")
        genericbodyforce1d = PhysicsSpace.wrapObjs["GenericBodyForce1d"]
        genericbodyforce2d = PhysicsSpace.wrapObjs["GenericBodyForce2d"]
        genericbodyforce3d = PhysicsSpace.wrapObjs["GenericBodyForce3d"]

        # Expose types.
        self.NBodyGravity1d = addObject(space, "NBodyGravity1d", allow_subclassing=True, parent=genericbodyforce1d)
        self.NBodyGravity2d = addObject(space, "NBodyGravity2d", allow_subclassing=True, parent=genericbodyforce2d)
        self.NBodyGravity3d = addObject(space, "NBodyGravity3d", allow_subclassing=True, parent=genericbodyforce3d)

        self.OctTreeGravity = addObject(space, "OctTreeGravity", allow_subclassing=True, parent=genericbodyforce3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings for all objects.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.generateNBodyGravityBindings(self.NBodyGravity1d, 1)
        self.generateNBodyGravityBindings(self.NBodyGravity2d, 2)
        self.generateNBodyGravityBindings(self.NBodyGravity3d, 3)

        self.generateOctTreeGravityBindings(self.OctTreeGravity)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["GravitySpace"]

    #---------------------------------------------------------------------------
    # Bindings (NBodyGravity)
    #---------------------------------------------------------------------------
    def generateNBodyGravityBindings(self, x, ndim):

        # Object names.
        me = "Spheral::GravitySpace::NBodyGravity%id" % ndim
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

        # Constructors.
        x.add_constructor([param("double", "plummerSofteningLength"),
                           param("double", "maxDeltaVelocity"),
                           param("double", "G")])

        # Methods.
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                                  constrefparam(state, "state"),
                                                  constrefparam(derivatives, "derivatives"),
                                                  param("double", "time")],
                     is_const=True, is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("softeningLength", "double", getter="softeningLength", setter="softeningLength")

        return

    #---------------------------------------------------------------------------
    # Bindings (OctTreeGravity)
    #---------------------------------------------------------------------------
    def generateOctTreeGravityBindings(self, x):

        # Object names.
        ndim = 3
        me = "Spheral::GravitySpace::OctTreeGravity"
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

        # Constructors.
        x.add_constructor([param("double", "G"),
                           param("double", "opening"),
                           param("double", "softeningLength"),
                           param("double", "ftimestep")])

        # Methods.
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                                  constrefparam(state, "state"),
                                                  constrefparam(derivatives, "derivatives"),
                                                  param("double", "time")],
                     is_const=True, is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("initialize", None, [param("const double", "time"),
                                          param("const double", "dt"),
                                          constrefparam(database, "dataBase"),
                                          refparam(state, "state"),
                                          refparam(derivatives, "derivs")], is_virtual=True)
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("opening", "double", getter="opening", setter="opening")
        x.add_instance_attribute("softeningLength", "double", getter="softeningLength", setter="softeningLength")
        x.add_instance_attribute("ftimestep", "double", getter="ftimestep", setter="ftimestep")
        x.add_instance_attribute("xmin", vector, getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", vector, getter="xmax", is_const=True)

        return
