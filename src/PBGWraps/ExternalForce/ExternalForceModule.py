from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ExternalForce:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/ExternalForceTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        for dim in self.dims:
            exec('''
genericbodyforce%(dim)id = findObject(space, "GenericBodyForce%(dim)id")

# Expose types.
self.PointPotential%(dim)id = addObject(space, "PointPotential%(dim)id", allow_subclassing=True, parent=genericbodyforce%(dim)id)
self.ConstantAcceleration%(dim)id = addObject(space, "ConstantAcceleration%(dim)id", allow_subclassing=True, parent=genericbodyforce%(dim)id)
self.LinearAcceleration%(dim)id = addObject(space, "LinearAcceleration%(dim)id", allow_subclassing=True, parent=genericbodyforce%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.addPointPotentialMethods(self.PointPotential%(dim)id, %(dim)i)
self.addConstantAccelerationMethods(self.ConstantAcceleration%(dim)id, %(dim)i)
self.addLinearAccelerationMethods(self.LinearAcceleration%(dim)id, %(dim)i)
''' % {"dim" : dim})

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
        me = "Spheral::ConstantAcceleration%id" % ndim
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
        const_ref_return_value(x, me, "%s::flags" % me, intfield, [], "flags")

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

        # Constructors.
        x.add_constructor([param("double", "a0"),
                           param("double", "aslope")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Attributes.
        x.add_instance_attribute("a0", "double", getter="a0", is_const=True)
        x.add_instance_attribute("aslope", "double", getter="aslope", is_const=True)

        return
