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
        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        space = Spheral.add_cpp_namespace("GravitySpace")

        for dim in self.dims:
            exec('''
genericbodyforce%(dim)id = findObject(PhysicsSpace, "GenericBodyForce%(dim)id")

# Expose types.
self.NBodyGravity%(dim)id = addObject(space, "NBodyGravity%(dim)id", allow_subclassing=True, parent=genericbodyforce%(dim)id)
''' % {"dim" : dim})

        if 2 in self.dims:
            self.QuadTreeGravity = addObject(space, "QuadTreeGravity", allow_subclassing=True, parent=genericbodyforce2d)
        if 3 in self.dims:
            self.OctTreeGravity =  addObject(space, "OctTreeGravity", allow_subclassing=True, parent=genericbodyforce3d)

        self.GravityTimeStepType = space.add_enum("GravityTimeStepType", ["GravityTimeStepType::AccelerationRatio", "GravityTimeStepType::DynamicalTime"])

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
            self.generateTreeGravityBindings(self.QuadTreeGravity, "Spheral::GravitySpace::QuadTreeGravity", 2)
        if 3 in self.dims:
            self.generateTreeGravityBindings(self.OctTreeGravity, "Spheral::GravitySpace::OctTreeGravity", 3)

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

        # Wrap the generic physics methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("softeningLength", "double", getter="softeningLength", setter="softeningLength")

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
                           param("double", "softeningLength"),
                           param("double", "opening", default_value="0.5"),
                           param("double", "ftimestep", default_value="0.1"),
                           param("GravityTimeStepType", "timeStepChoice", default_value="Spheral::GravitySpace::GravityTimeStepType::AccelerationRatio")])

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
