from pybindgen import *

from PBGutils import *
from ref_return_value import *

from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Hydro:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/HydroTypes.hh"' % srcdir)
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("PhysicsSpace")

        # Expose types.
        self.HydroFieldNames = addObject(Spheral, "HydroFieldNames")

        for dim in self.dims:
            exec('''
physics%(dim)id = findObject(space, "Physics%(dim)id")
self.ThirdMomentHourglassControl%(dim)id = addObject(space, "ThirdMomentHourglassControl%(dim)id", allow_subclassing=True, parent=physics%(dim)id)
''' % {"dim" : dim})

        if 1 in self.dims:
            self.SecondMomentHourglassControl1d = addObject(space, "SecondMomentHourglassControl1d", allow_subclassing=True, parent=physics1d)
        if 2 in self.dims:
            self.SecondMomentHourglassControl2d = addObject(space, "SecondMomentHourglassControl2d", allow_subclassing=True, parent=physics2d)

        # self.VoronoiHourglassControl1d = addObject(space, "VoronoiHourglassControl1d", allow_subclassing=True, parent=physics1d)
        # self.VoronoiHourglassControl2d = addObject(space, "VoronoiHourglassControl2d", allow_subclassing=True, parent=physics2d)
        # self.VoronoiHourglassControl3d = addObject(space, "VoronoiHourglassControl3d", allow_subclassing=True, parent=physics3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateHydroFieldNamesBindings(self.HydroFieldNames)

        if 1 in self.dims:
            self.generateSecondMomentHourglassControlBindings(self.SecondMomentHourglassControl1d, 1)
        if 2 in self.dims:
            self.generateSecondMomentHourglassControlBindings(self.SecondMomentHourglassControl2d, 2)

        for dim in self.dims:
            exec('''
self.generateThirdMomentHourglassControlBindings(self.ThirdMomentHourglassControl%(dim)id, %(dim)i)
''' % {"dim" : dim})

        # self.generateVoronoiHourglassControlBindings(self.VoronoiHourglassControl1d, 1)
        # self.generateVoronoiHourglassControlBindings(self.VoronoiHourglassControl2d, 2)
        # self.generateVoronoiHourglassControlBindings(self.VoronoiHourglassControl3d, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # Bindings (HydroFieldNames).
    #---------------------------------------------------------------------------
    def generateHydroFieldNamesBindings(self, x):
        x.add_static_attribute("mass", "std::string",  is_const=True)
        x.add_static_attribute("position", "std::string",  is_const=True)
        x.add_static_attribute("velocity", "std::string",  is_const=True)
        x.add_static_attribute("H", "std::string",  is_const=True)
        x.add_static_attribute("work", "std::string",  is_const=True)
        x.add_static_attribute("velocityGradient", "std::string",  is_const=True)
        x.add_static_attribute("internalVelocityGradient", "std::string",  is_const=True)
        x.add_static_attribute("massDensity", "std::string",  is_const=True)
        x.add_static_attribute("normalization", "std::string",  is_const=True)
        x.add_static_attribute("specificThermalEnergy", "std::string",  is_const=True)
        x.add_static_attribute("maxViscousPressure", "std::string",  is_const=True)
        x.add_static_attribute("effectiveViscousPressure", "std::string",  is_const=True)
        x.add_static_attribute("viscousWork", "std::string",  is_const=True)
        x.add_static_attribute("XSPHDeltaV", "std::string",  is_const=True)
        x.add_static_attribute("XSPHWeightSum", "std::string",  is_const=True)
        x.add_static_attribute("Hsmooth", "std::string",  is_const=True)
        x.add_static_attribute("massFirstMoment", "std::string",  is_const=True)
        x.add_static_attribute("massSecondMoment", "std::string",  is_const=True)
        x.add_static_attribute("weightedNeighborSum", "std::string",  is_const=True)
        x.add_static_attribute("pressure", "std::string",  is_const=True)
        x.add_static_attribute("temperature", "std::string",  is_const=True)
        x.add_static_attribute("soundSpeed", "std::string",  is_const=True)
        x.add_static_attribute("pairAccelerations", "std::string",  is_const=True)
        x.add_static_attribute("pairWork", "std::string",  is_const=True)
        x.add_static_attribute("omegaGradh", "std::string",  is_const=True)
        x.add_static_attribute("numberDensitySum", "std::string",  is_const=True)
        x.add_static_attribute("timeStepMask", "std::string",  is_const=True)
        x.add_static_attribute("m0_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("m1_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("m2_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("A0_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("A_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("B_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("gradA_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("gradB_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("M_CRKSPH", "std::string",  is_const=True)
        x.add_static_attribute("volume", "std::string",  is_const=True)
        x.add_static_attribute("linearMomentum", "std::string",  is_const=True)
        x.add_static_attribute("totalEnergy", "std::string",  is_const=True)
        x.add_static_attribute("mesh", "std::string",  is_const=True)
        x.add_static_attribute("hourglassMask", "std::string",  is_const=True)
        x.add_static_attribute("faceVelocity", "std::string",  is_const=True)
        x.add_static_attribute("faceForce", "std::string",  is_const=True)
        x.add_static_attribute("faceMass", "std::string",  is_const=True)
        x.add_static_attribute("polyvols", "std::string",  is_const=True)
        x.add_static_attribute("massDensityGradient", "std::string",  is_const=True)
        return

    #---------------------------------------------------------------------------
    # Bindings (SecondMomentHourglassControl).
    #---------------------------------------------------------------------------
    def generateSecondMomentHourglassControlBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::SecondMomentHourglassControl%id" % ndim
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

        # Constructors.
        x.add_constructor([constrefparam(tablekernel, "W"),
                           param("double", "multiplier", default_value="0.05"),
                           param("double", "maxAccelerationFactor", default_value="0.001")])

        # Attributes.
        x.add_instance_attribute("multiplier", "double", getter="multiplier", setter="multiplier")
        x.add_instance_attribute("maxAccelerationFactor", "double", getter="maxAccelerationFactor", setter="maxAccelerationFactor")
        x.add_instance_attribute("acceleration", vectorfieldlist, getter="acceleration", is_const=True)

        # Override the abstract interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Bindings (ThirdMomentHourglassControl).
    #---------------------------------------------------------------------------
    def generateThirdMomentHourglassControlBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::ThirdMomentHourglassControl%id" % ndim
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

        # Constructors.
        x.add_constructor([constrefparam(database, "dataBase"),
                           constrefparam(tablekernel, "W"),
                           param("double", "multiplier", default_value="0.05"),
                           param("double", "maxAccelerationFactor", default_value="0.01")])

        # Attributes.
        x.add_instance_attribute("multiplier", "double", getter="multiplier", setter="multiplier")
        x.add_instance_attribute("maxAccelerationFactor", "double", getter="maxAccelerationFactor", setter="maxAccelerationFactor")
        x.add_instance_attribute("thirdMoment", thirdranktensorfieldlist, getter="thirdMoment", is_const=True)

        # Override the abstract interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Bindings (VoronoiHourglassControl).
    #---------------------------------------------------------------------------
    def generateVoronoiHourglassControlBindings(self, x, ndim):

        # Object names.
        me = "Spheral::PhysicsSpace::VoronoiHourglassControl%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(tablekernel, "W"),
                           param("unsigned int", "order", default_value="1"),
                           param("unsigned int", "limiter", default_value="1"),
                           param("double", "fraction", default_value="0.5"),
                           param(intfieldlist, "mask", default_value=("%s(FieldSpace::Copy)" % intfieldlist))])

        # Methods.
        const_ref_return_value(x, me, "%s::mask" % me, intfieldlist, [], "mask")
        const_ref_return_value(x, me, "%s::kernel" % me, tablekernel, [], "kernel")
        const_ref_return_value(x, me, "%s::gradRho" % me, vectorfieldlist, [], "gradRho")
        const_ref_return_value(x, me, "%s::A" % me, scalarfieldlist, [], "A")
        const_ref_return_value(x, me, "%s::B" % me, vectorfieldlist, [], "B")
        const_ref_return_value(x, me, "%s::C" % me, vectorfieldlist, [], "C")
        const_ref_return_value(x, me, "%s::D" % me, tensorfieldlist, [], "D")
        const_ref_return_value(x, me, "%s::gradA" % me, vectorfieldlist, [], "gradA")
        const_ref_return_value(x, me, "%s::gradB" % me, tensorfieldlist, [], "gradB")
        const_ref_return_value(x, me, "%s::weight" % me, scalarfieldlist, [], "weight")

        # Attributes.
        x.add_instance_attribute("order", "unsigned int", getter="order", setter="order")
        x.add_instance_attribute("limiter", "unsigned int", getter="limiter", setter="limiter")
        x.add_instance_attribute("fraction", "double", getter="fraction", setter="fraction")

        # Override the abstract interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        return
