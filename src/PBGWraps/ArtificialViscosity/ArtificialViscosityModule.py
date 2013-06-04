from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ArtificialViscosity:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"ArtificialViscosity/ArtificialViscosityTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("ArtificialViscositySpace")

        # Expose types.
        self.ArtificialViscosity1d = addObject(space, "ArtificialViscosity1d", allow_subclassing=True)
        self.ArtificialViscosity2d = addObject(space, "ArtificialViscosity2d", allow_subclassing=True)
        self.ArtificialViscosity3d = addObject(space, "ArtificialViscosity3d", allow_subclassing=True)

        self.MonaghanGingoldViscosity1d = addObject(space, "MonaghanGingoldViscosity1d", allow_subclassing=True, parent=self.ArtificialViscosity1d)
        self.MonaghanGingoldViscosity2d = addObject(space, "MonaghanGingoldViscosity2d", allow_subclassing=True, parent=self.ArtificialViscosity2d)
        self.MonaghanGingoldViscosity3d = addObject(space, "MonaghanGingoldViscosity3d", allow_subclassing=True, parent=self.ArtificialViscosity3d)

        self.TensorMonaghanGingoldViscosity1d = addObject(space, "TensorMonaghanGingoldViscosity1d", allow_subclassing=True, parent=self.ArtificialViscosity1d)
        self.TensorMonaghanGingoldViscosity2d = addObject(space, "TensorMonaghanGingoldViscosity2d", allow_subclassing=True, parent=self.ArtificialViscosity2d)
        self.TensorMonaghanGingoldViscosity3d = addObject(space, "TensorMonaghanGingoldViscosity3d", allow_subclassing=True, parent=self.ArtificialViscosity3d)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.addArtificialViscosityMethods(self.ArtificialViscosity1d, 1)
        self.addArtificialViscosityMethods(self.ArtificialViscosity2d, 2)
        self.addArtificialViscosityMethods(self.ArtificialViscosity3d, 3)

        self.addMonaghanGingoldViscosityMethods(self.MonaghanGingoldViscosity1d, 1)
        self.addMonaghanGingoldViscosityMethods(self.MonaghanGingoldViscosity2d, 2)
        self.addMonaghanGingoldViscosityMethods(self.MonaghanGingoldViscosity3d, 3)

        self.addTensorMonaghanGingoldViscosityMethods(self.TensorMonaghanGingoldViscosity1d, 1)
        self.addTensorMonaghanGingoldViscosityMethods(self.TensorMonaghanGingoldViscosity2d, 2)
        self.addTensorMonaghanGingoldViscosityMethods(self.TensorMonaghanGingoldViscosity3d, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["ArtificialViscositySpace"]

    #---------------------------------------------------------------------------
    # Add methods to the ArtificialViscosity.
    #---------------------------------------------------------------------------
    def addArtificialViscosityMethods(self, x, ndim):

        me = "ArtificialViscosity%id" % ndim

        # External objects.
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0")])

        # Attributes.
        x.add_instance_attribute("Cl", "double", getter="Cl", setter="Cl")
        x.add_instance_attribute("Cq", "double", getter="Cq", setter="Cq")
        x.add_instance_attribute("balsaraShearCorrection", "bool", getter="balsaraShearCorrection", setter="balsaraShearCorrection")
        x.add_instance_attribute("limiter", "bool", getter="limiter", setter="limiter")
        x.add_instance_attribute("epsilon2", "double", getter="epsilon2", setter="epsilon2")
        x.add_instance_attribute("negligibleSoundSpeed", "double", getter="negligibleSoundSpeed", setter="negligibleSoundSpeed")
        x.add_instance_attribute("csMultiplier", "double", getter="csMultiplier", setter="csMultiplier")
        x.add_instance_attribute("energyMultiplier", "double", getter="energyMultiplier", setter="energyMultiplier")

        # Virtual methods for all Qs.
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)

        # Methods.
        x.add_method("curlVelocityMagnitude", "double", [refparam(tensor, "DvDx")], is_const=True)
        x.add_method("shearMultiplier", scalarfieldlist, [], is_const=True)
        x.add_method("sigma", tensorfieldlist, [], is_const=True)
        x.add_method("gradDivVelocity", vectorfieldlist, [], is_const=True)
        x.add_method("calculateLimiter", tensor, [refparam(vector, "vi"),
                                                  refparam(vector, "vj"),
                                                  param("double", "ci"),
                                                  param("double", "cj"),
                                                  param("double", "hi"),
                                                  param("double", "hj"),
                                                  param("int", "nodeListID"),
                                                  param("int", "nodeID")],
                     is_const = True)
        x.add_method("shockDirection", vector, [param("double", "ci"),
                                                param("double", "hi"),
                                                param("int", "nodeListID"),
                                                param("int", "nodeID")],
                     is_const = True)
        x.add_method("sigmaij", tensor, [refparam(vector, "rji"),
                                         refparam(vector, "rjiUnit"),
                                         refparam(vector, "vji"),
                                         refparam("double", "hi2"),
                                         param("int", "nodeListID"),
                                         param("int", "nodeID")],
                     is_const = True)
        x.add_method("calculateSigma", "bool", [], is_const=True, visibility="protected")
        x.add_method("calculateSigma", None, [param("bool", "value")], visibility="protected")

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # Add methods to the MonaghanGingoldViscosity.
    #---------------------------------------------------------------------------
    def addMonaghanGingoldViscosityMethods(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Add methods to the TensorMonaghanGingoldViscosity.
    #---------------------------------------------------------------------------
    def addTensorMonaghanGingoldViscosityMethods(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Add pure virtual methods to the ArtificialViscosity.
    #---------------------------------------------------------------------------
    def addArtificialViscosityVirtualMethods(self, x, ndim, pureVirtual):

        # External objects.
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        pair_tensor_tensor = "pair_Tensor%id_Tensor%id" % (ndim, ndim)

        # Abstract Q interface.
        x.add_method("Piij", pair_tensor_tensor, [param("unsigned int", "nodeListi"),
                                                  param("unsigned int", "i"),
                                                  param("unsigned int", "nodeListj"),
                                                  param("unsigned int", "j"),
                                                  constrefparam(vector, "xi"),
                                                  constrefparam(vector, "etai"),
                                                  constrefparam(vector, "vi"),
                                                  param("double", "rhoi"),
                                                  param("double", "csi"),
                                                  constrefparam(symtensor, "Hi"),
                                                  constrefparam(vector, "xj"),
                                                  constrefparam(vector, "etaj"),
                                                  constrefparam(vector, "vj"),
                                                  param("double", "rhoj"),
                                                  param("double", "csj"),
                                                  constrefparam(symtensor, "Hj")],
                     is_const = True,
                     is_virtual = True,
                     is_pure_virtual = pureVirtual)

        return

