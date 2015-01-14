from pybindgen import *

import sys
srcdir = sys.argv[2]

from PBGutils import *
from ref_return_value import *

sys.path.append("%s/Physics" % srcdir)
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ArtificialViscosity:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"%s/ArtificialViscosity/ArtificialViscosityTypes.hh"' % srcdir)
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("ArtificialViscositySpace")
        physics = Spheral.add_cpp_namespace("PhysicsSpace")
        
        Physics1d = findObject(physics,"Physics1d")
        Physics2d = findObject(physics,"Physics2d")
        Physics3d = findObject(physics,"Physics3d")

        # Expose types.
        self.dimSet = (1, 2, 3)
        for dim in self.dimSet:
            exec('''
self.ArtificialViscosity%(dim)id = addObject(space, "ArtificialViscosity%(dim)id", allow_subclassing=True)
self.MonaghanGingoldViscosity%(dim)id = addObject(space, "MonaghanGingoldViscosity%(dim)id", allow_subclassing=True, parent=self.ArtificialViscosity%(dim)id)
self.CSPHMonaghanGingoldViscosity%(dim)id = addObject(space, "CSPHMonaghanGingoldViscosity%(dim)id", allow_subclassing=True, parent=self.MonaghanGingoldViscosity%(dim)id)
self.MorrisMonaghanReducingViscosity%(dim)id = addObject(space, "MorrisMonaghanReducingViscosity%(dim)id", allow_subclassing=True, parent=Physics%(dim)id)
self.TensorMonaghanGingoldViscosity%(dim)id = addObject(space, "TensorMonaghanGingoldViscosity%(dim)id", allow_subclassing=True, parent=self.ArtificialViscosity%(dim)id)
self.FiniteVolumeViscosity%(dim)id = addObject(space, "FiniteVolumeViscosity%(dim)id", allow_subclassing=True, parent=self.ArtificialViscosity%(dim)id)
self.TensorSVPHViscosity%(dim)id = addObject(space, "TensorSVPHViscosity%(dim)id", allow_subclassing=True, parent=self.ArtificialViscosity%(dim)id)
self.TensorCSPHViscosity%(dim)id = addObject(space, "TensorCSPHViscosity%(dim)id", allow_subclassing=True, parent=self.TensorMonaghanGingoldViscosity%(dim)id)
self.VonNeumanViscosity%(dim)id = addObject(space, "VonNeumanViscosity%(dim)id", allow_subclassing=True, parent=self.ArtificialViscosity%(dim)id)
''' % {"dim" : dim})
        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        for dim in self.dimSet:
            exec('''
self.addArtificialViscosityMethods(self.ArtificialViscosity%(dim)id, %(dim)i)
self.addMonaghanGingoldViscosityMethods(self.MonaghanGingoldViscosity%(dim)id, %(dim)i)
self.addCSPHMonaghanGingoldViscosityMethods(self.CSPHMonaghanGingoldViscosity%(dim)id, %(dim)i)
self.addMorrisMonaghanReducingViscosityMethods(self.MorrisMonaghanReducingViscosity%(dim)id, %(dim)i)
self.addTensorMonaghanGingoldViscosityMethods(self.TensorMonaghanGingoldViscosity%(dim)id, %(dim)i)
self.addFiniteVolumeViscosityMethods(self.FiniteVolumeViscosity%(dim)id, %(dim)i)
self.addTensorSVPHViscosityMethods(self.TensorSVPHViscosity%(dim)id, %(dim)i)
self.addTensorCSPHViscosityMethods(self.TensorCSPHViscosity%(dim)id, %(dim)i)
self.addVonNeumanViscosityMethods(self.VonNeumanViscosity%(dim)id, %(dim)i)
''' % {"dim" : dim})
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

        me = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim

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
        x.add_instance_attribute("reducingViscosityCorrection", "bool", getter="reducingViscosityCorrection", setter="reducingViscosityCorrection")
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
        
        x.add_method("reducingViscosityMultiplierQ", scalarfieldlist, [])
        x.add_method("reducingViscosityMultiplierL", scalarfieldlist, [])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # Add methods to the MonaghanGingoldViscosity.
    #---------------------------------------------------------------------------
    def addMonaghanGingoldViscosityMethods(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0"),
                           param("bool", "linearInExpansion", default_value="false"),
                           param("bool", "quadraticInExpansion", default_value="false")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        # Attributes
        x.add_instance_attribute("linearInExpansion", "bool", getter="linearInExpansion", setter="linearInExpansion")
        x.add_instance_attribute("quadraticInExpansion", "bool", getter="quadraticInExpansion", setter="quadraticInExpansion")

        return
    
    #---------------------------------------------------------------------------
    # Add methods to the CSPHMonaghanGingoldViscosity.
    #---------------------------------------------------------------------------
    def addCSPHMonaghanGingoldViscosityMethods(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0"),
                           param("bool", "linearInExpansion", default_value="false"),
                           param("bool", "quadraticInExpansion", default_value="false")])

        # Add the local methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        return
    
    #---------------------------------------------------------------------------
    # Add methods to the MorrsMonaghanReducingViscosity.
    #---------------------------------------------------------------------------
    def addMorrisMonaghanReducingViscosityMethods(self, x, ndim):
        
        me = "Spheral::ArtificialViscositySpace::MorrisMonaghanReducingViscosity%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        
        # Constructors.
        x.add_constructor([refparam(artificialviscosity,"q"),
                           param("double", "nhQ", default_value="5.0"),
                           param("double", "nhL", default_value="10.0"),
                           param("double", "aMin", default_value="0.1"),
                           param("double", "aMax", default_value="2.0")])
                           
        # Add the abstract methods.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Attributes
        x.add_instance_attribute("nhQ", "double", getter="nhQ", setter="nhQ")
        x.add_instance_attribute("nhL", "double", getter="nhL", setter="nhL")
        x.add_instance_attribute("aMin", "double", getter="aMin", setter="aMin")
        x.add_instance_attribute("aMax", "double", getter="aMax", setter="aMax")
        
        # Methods.
        const_ref_return_value(x, me, "%s::DrvAlphaDtQ" % me, scalarfieldlist, [], "DrvAlphaDtQ")
        const_ref_return_value(x, me, "%s::DrvAlphaDtL" % me, scalarfieldlist, [], "DrvAlphaDtL")
    
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
    # Add methods to the FiniteVolumeViscosity.
    #---------------------------------------------------------------------------
    def addFiniteVolumeViscosityMethods(self, x, ndim):

        me = "Spheral::ArtificialViscositySpace::FiniteVolumeViscosity%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0"),
                           param("bool", "scalar", default_value="false")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        # Attributes
        x.add_instance_attribute("scalar", "bool", getter="scalar", is_const=True)
        const_ref_return_value(x, me, "%s::DvDx" % me, tensorfieldlist, [], "DvDx")

        return

    #---------------------------------------------------------------------------
    # Add methods to the TensorSVPHViscosity.
    #---------------------------------------------------------------------------
    def addTensorSVPHViscosityMethods(self, x, ndim):

        me = "Spheral::ArtificialViscositySpace::TensorSVPHViscosity%id" % ndim
        vector_of_tensor = "vector_of_Tensor%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0"),
                           param("double", "fslice", default_value="0.5")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        # Attributes.
        x.add_instance_attribute("fslice", "double", getter="fslice", setter="fslice")
        const_ref_return_value(x, me, "%s::DvDx" % me, vector_of_tensor, [], "DvDx")
        const_ref_return_value(x, me, "%s::shearCorrection" % me, "vector_of_double", [], "shearCorrection")
        const_ref_return_value(x, me, "%s::Qface" % me, vector_of_tensor, [], "Qface")

        return

    #---------------------------------------------------------------------------
    # Add methods to the TensorCSPHViscosity.
    #---------------------------------------------------------------------------
    def addTensorCSPHViscosityMethods(self, x, ndim):

        me = "Spheral::ArtificialViscositySpace::TensorCSPHViscosity%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Add methods to the VonNeumanViscosity.
    #---------------------------------------------------------------------------
    def addVonNeumanViscosityMethods(self, x, ndim):

        me = "Spheral::ArtificialViscositySpace::VonNeumanViscosity%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim

        # Constructors.
        x.add_constructor([param("double", "Clinear", default_value="1.0"),
                           param("double", "Cquadratic", default_value="1.0")])

        # Add the abstract methods.
        self.addArtificialViscosityVirtualMethods(x, ndim, False)

        # Attributes
        const_ref_return_value(x, me, "%s::viscousEnergy" % me, scalarfieldlist, [], "DvDx")

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

