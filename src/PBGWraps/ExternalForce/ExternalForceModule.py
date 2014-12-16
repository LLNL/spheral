from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

sys.path.append("../Physics")
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ExternalForce:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"ExternalForce/ExternalForceTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("PhysicsSpace")

        genericbodyforce1d = findObject(space, "GenericBodyForce1d")
        genericbodyforce2d = findObject(space, "GenericBodyForce2d")
        genericbodyforce3d = findObject(space, "GenericBodyForce3d")

        # Expose types.
        self.PointPotential1d = addObject(space, "PointPotential1d", allow_subclassing=True, parent=genericbodyforce1d)
        self.PointPotential2d = addObject(space, "PointPotential2d", allow_subclassing=True, parent=genericbodyforce2d)
        self.PointPotential3d = addObject(space, "PointPotential3d", allow_subclassing=True, parent=genericbodyforce3d)

        self.ConstantAcceleration1d = addObject(space, "ConstantAcceleration1d", allow_subclassing=True, parent=genericbodyforce1d)
        self.ConstantAcceleration2d = addObject(space, "ConstantAcceleration2d", allow_subclassing=True, parent=genericbodyforce2d)
        self.ConstantAcceleration3d = addObject(space, "ConstantAcceleration3d", allow_subclassing=True, parent=genericbodyforce3d)

        self.LinearAcceleration1d = addObject(space, "LinearAcceleration1d", allow_subclassing=True, parent=genericbodyforce1d)
        self.LinearAcceleration2d = addObject(space, "LinearAcceleration2d", allow_subclassing=True, parent=genericbodyforce2d)
        self.LinearAcceleration3d = addObject(space, "LinearAcceleration3d", allow_subclassing=True, parent=genericbodyforce3d)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addPointPotentialMethods(self.PointPotential1d, 1)
        self.addPointPotentialMethods(self.PointPotential2d, 2)
        self.addPointPotentialMethods(self.PointPotential3d, 3)

        self.addConstantAccelerationMethods(self.ConstantAcceleration1d, 1)
        self.addConstantAccelerationMethods(self.ConstantAcceleration2d, 2)
        self.addConstantAccelerationMethods(self.ConstantAcceleration3d, 3)

        self.addLinearAccelerationMethods(self.LinearAcceleration1d, 1)
        self.addLinearAccelerationMethods(self.LinearAcceleration2d, 2)
        self.addLinearAccelerationMethods(self.LinearAcceleration3d, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # Add PointPotential methods.
    #---------------------------------------------------------------------------
    def addPointPotentialMethods(self, x, ndim):

        # Object names.
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
                           param("double", "mass"),
                           param("double", "coreRadius"),
                           constrefparam(vector, "origin")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("specificPotential", "double", [constrefparam(vector, "r")], is_const=True)

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", setter="setG")
        x.add_instance_attribute("mass", "double", getter="mass", setter="setMass")
        x.add_instance_attribute("coreRadius", "double", getter="coreRadius", setter="setCoreRadius")
        x.add_instance_attribute("origin", vector, getter="origin", setter="setOrigin")
        x.add_instance_attribute("deltaPotentialFraction", "double", getter="deltaPotentialFraction", setter="setDeltaPotentialFraction")

        return

    #---------------------------------------------------------------------------
    # Add ConstantAcceleration methods.
    #---------------------------------------------------------------------------
    def addConstantAccelerationMethods(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::ConstantAcceleration%id" % ndim
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
        x.add_constructor([param(vector, "a0"),
                           constrefparam(nodelist, "nodeList"),
                           constrefparam("vector_of_int", "indices")])
        x.add_constructor([param(vector, "a0"),
                           constrefparam(nodelist, "nodeList")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
        const_ref_return_value(x, me, "%s::indices" % me, "vector_of_int", [], "indices")

        # Attributes.
        x.add_instance_attribute("a0", vector, getter="a0", is_const=True)

        return

    #---------------------------------------------------------------------------
    # Add LinearAcceleration methods.
    #---------------------------------------------------------------------------
    def addLinearAccelerationMethods(self, x, ndim):

        # Object names.
        me = "LinearAcceleration%id" % ndim
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
        x.add_constructor([param("double", "a0"),
                           param("double", "aslope")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Attributes.
        x.add_instance_attribute("a0", "double", getter="a0", is_const=True)
        x.add_instance_attribute("aslope", "double", getter="aslope", is_const=True)

        return
