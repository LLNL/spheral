from pybindgen import *

from ref_return_value import *
from PhysicsModule import generatePhysicsVirtualBindings
from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Strength:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/StrengthTypes.hh"' % srcdir)
    
        # Namespace.
        SolidSpheral = mod.add_cpp_namespace("Spheral")
        space = SolidSpheral.add_cpp_namespace("SolidMaterial")

        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        NodeSpace = Spheral.add_cpp_namespace("NodeSpace")

        self.SolidFieldNames = addObject(SolidSpheral, "SolidFieldNames")

        for dim in self.dims:
            exec('''
Physics%(dim)id = findObject(PhysicsSpace, "Physics%(dim)id")
FluidNodeList%(dim)id = findObject(NodeSpace, "FluidNodeList%(dim)id")
self.SolidNodeList%(dim)id = addObject(space, "SolidNodeList%(dim)id", parent=FluidNodeList%(dim)id, allow_subclassing=True)
self.vector_of_SolidNodeList%(dim)id = addObject(mod, "vector_of_SolidNodeList%(dim)id", allow_subclassing=True)
self.vector_of_SolidNodeList%(dim)id_iterator = addObject(mod, "vector_of_SolidNodeList%(dim)id_iterator", allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateSolidFieldNamesBindings(self.SolidFieldNames)

        for dim in self.dims:
            exec('''
self.generateSolidNodeListBindings(self.SolidNodeList%(dim)id, %(dim)i)
generateStdVectorBindings(self.vector_of_SolidNodeList%(dim)id, "Spheral::SolidMaterial::SolidNodeList%(dim)id*", "vector_of_SolidNodeList%(dim)id")
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["SolidMaterial"]

    #---------------------------------------------------------------------------
    # Bindings (SolidFieldNames).
    #---------------------------------------------------------------------------
    def generateSolidFieldNamesBindings(self, x):
        x.add_static_attribute("deviatoricStress", "std::string",  is_const=True)
        x.add_static_attribute("plasticStrain", "std::string",  is_const=True)
        x.add_static_attribute("scalarDamage", "std::string",  is_const=True)
        x.add_static_attribute("tensorDamage", "std::string",  is_const=True)
        x.add_static_attribute("effectiveTensorDamage", "std::string",  is_const=True)
        x.add_static_attribute("damageGradient", "std::string",  is_const=True)
        x.add_static_attribute("damageHat", "std::string",  is_const=True)
        x.add_static_attribute("strain", "std::string",  is_const=True)
        x.add_static_attribute("strainTensor", "std::string",  is_const=True)
        x.add_static_attribute("effectiveStrainTensor", "std::string",  is_const=True)
        x.add_static_attribute("bulkModulus", "std::string",  is_const=True)
        x.add_static_attribute("shearModulus", "std::string",  is_const=True)
        x.add_static_attribute("YoungsModulus", "std::string",  is_const=True)
        x.add_static_attribute("longitudinalSoundSpeed", "std::string",  is_const=True)
        x.add_static_attribute("yieldStrength", "std::string",  is_const=True)
        x.add_static_attribute("effectiveFlaws", "std::string",  is_const=True)
        x.add_static_attribute("fragmentIDs", "std::string",  is_const=True)
        return

    #---------------------------------------------------------------------------
    # SolidNodeList
    #---------------------------------------------------------------------------
    def generateSolidNodeListBindings(self, x, ndim):

        me = "Spheral::SolidMaterial::SolidNodeList%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::Material::EquationOfState%id" % ndim
        strengthmodel = "Spheral::SolidMaterial::StrengthModel%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([param("std::string", "name"),
                           refparam(equationofstate, "eos"),
                           refparam(strengthmodel, "strength"),
                           param("int", "numInternal", default_value="0"),
                           param("int", "numGhost", default_value="0"),
                           param("double", "hmin", default_value="0.0"),
                           param("double", "hmax", default_value="1.0e100"),
                           param("double", "hminratio", default_value="0.1"),
                           param("double", "nPerh", default_value="2.01"),
                           param("int", "maxNumNeighbors", default_value="500"),
                           param("double", "rhoMin", default_value="1.0e-10"),
                           param("double", "rhoMax", default_value="1.0e100")])

        # Methods.
        x.add_method("soundSpeed", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("bulkModulus", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("shearModulus", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("yieldStrength", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)

        const_ref_return_value(x, me, "%s::deviatoricStress" % me, symtensorfield, [], "deviatoricStress")
        const_ref_return_value(x, me, "%s::plasticStrain" % me, scalarfield, [], "plasticStrain")
        const_ref_return_value(x, me, "%s::plasticStrainRate" % me, scalarfield, [], "plasticStrain")
        const_ref_return_value(x, me, "%s::damage" % me, symtensorfield, [], "damage")
        const_ref_return_value(x, me, "%s::effectiveDamage" % me, symtensorfield, [], "effectiveDamage")
        const_ref_return_value(x, me, "%s::damageGradient" % me, vectorfield, [], "damageGradient")
        const_ref_return_value(x, me, "%s::fragmentIDs" % me, intfield, [], "fragmentIDs")
        const_ref_return_value(x, me, "%s::strengthModel" % me, strengthmodel, [], "strengthModel")

        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        return

