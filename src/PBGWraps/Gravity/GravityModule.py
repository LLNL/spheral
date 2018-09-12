from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Gravity:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/GravityTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        for dim in self.dims:
            exec('''
genericbodyforce%(dim)id = findObject(space, "GenericBodyForce%(dim)id")

# Expose types.
self.NBodyGravity%(dim)id = addObject(space, "NBodyGravity%(dim)id", allow_subclassing=True, parent=genericbodyforce%(dim)id)
''' % {"dim" : dim})

        if 2 in self.dims:
            self.QuadTreeGravity = addObject(space, "QuadTreeGravity", allow_subclassing=True, parent=genericbodyforce2d)
        if 3 in self.dims:
            self.OctTreeGravity =  addObject(space, "OctTreeGravity", allow_subclassing=True, parent=genericbodyforce3d)

        self.GravityTimeStepType = space.add_enum("GravityTimeStepType", [("AccelerationRatio", "Spheral::GravityTimeStepType::AccelerationRatio"),
                                                                          ("DynamicalTime", "Spheral::GravityTimeStepType::DynamicalTime")])

        return

    #---------------------------------------------------------------------------
    # Generate bindings for all objects.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        for dim in self.dims:
            exec('''
self.generateNBodyGravityBindings(self.NBodyGravity%(dim)id, %(dim)i)
''' % {"dim" : dim})

        if 2 in self.dims:
            self.generateTreeGravityBindings(self.QuadTreeGravity, "Spheral::QuadTreeGravity", 2)
        if 3 in self.dims:
            self.generateTreeGravityBindings(self.OctTreeGravity, "Spheral::OctTreeGravity", 3)

        return

    #---------------------------------------------------------------------------
    # Bindings (NBodyGravity)
    #---------------------------------------------------------------------------
    def generateNBodyGravityBindings(self, x, ndim):

        # Object names.
        me = "Spheral::NBodyGravity%id" % ndim
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
        x.add_constructor([param("double", "plummerSofteningLength"),
                           param("double", "maxDeltaVelocity"),
                           param("double", "G"),
                           param("bool", "compatibleVelocityUpdate", default_value="false")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("softeningLength", "double", getter="softeningLength", setter="softeningLength")
        x.add_instance_attribute("compatibleVelocityUpdate", "bool", getter="compatibleVelocityUpdate", setter="compatibleVelocityUpdate")

        return

    #---------------------------------------------------------------------------
    # Bindings (TreeGravity)
    #---------------------------------------------------------------------------
    def generateTreeGravityBindings(self, x, me, ndim):

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
                           param("double", "softeningLength"),
                           param("double", "opening", default_value="0.5"),
                           param("double", "ftimestep", default_value="0.1"),
                           param("GravityTimeStepType", "timeStepChoice", default_value="Spheral::GravityTimeStepType::AccelerationRatio")])

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("dumpTree", "std::string", [param("bool", "globalTree")], is_const=True)
        x.add_method("dumpTreeStatistics", "std::string", [param("bool", "globalTree")], is_const=True)
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("opening", "double", getter="opening", setter="opening")
        x.add_instance_attribute("softeningLength", "double", getter="softeningLength", setter="softeningLength")
        x.add_instance_attribute("ftimestep", "double", getter="ftimestep", setter="ftimestep")
        x.add_instance_attribute("timeStepChoice", "GravityTimeStepType", getter="timeStepChoice", setter="timeStepChoice")
        x.add_instance_attribute("xmin", vector, getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", vector, getter="xmax", is_const=True)

        return
