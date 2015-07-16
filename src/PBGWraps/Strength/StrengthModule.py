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
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/StrengthTypes.hh"' % srcdir)
    
        # Namespace.
        SolidSpheral = mod.add_cpp_namespace("Spheral")
        space = SolidSpheral.add_cpp_namespace("SolidMaterial")

        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        NodeSpace = Spheral.add_cpp_namespace("NodeSpace")

        Physics1d = findObject(PhysicsSpace, "Physics1d")
        Physics2d = findObject(PhysicsSpace, "Physics2d")
        Physics3d = findObject(PhysicsSpace, "Physics3d")

        FluidNodeList1d = findObject(NodeSpace, "FluidNodeList1d")
        FluidNodeList2d = findObject(NodeSpace, "FluidNodeList2d")
        FluidNodeList3d = findObject(NodeSpace, "FluidNodeList3d")

        self.SolidFieldNames = addObject(SolidSpheral, "SolidFieldNames")

        self.SolidNodeList1d = addObject(space, "SolidNodeList1d", parent=FluidNodeList1d, allow_subclassing=True)
        self.SolidNodeList2d = addObject(space, "SolidNodeList2d", parent=FluidNodeList2d, allow_subclassing=True)
        self.SolidNodeList3d = addObject(space, "SolidNodeList3d", parent=FluidNodeList3d, allow_subclassing=True)

        self.vector_of_SolidNodeList1d = addObject(mod, "vector_of_SolidNodeList1d", allow_subclassing=True)
        self.vector_of_SolidNodeList2d = addObject(mod, "vector_of_SolidNodeList2d", allow_subclassing=True)
        self.vector_of_SolidNodeList3d = addObject(mod, "vector_of_SolidNodeList3d", allow_subclassing=True)

        self.vector_of_SolidNodeList1d_iterator = addObject(mod, "vector_of_SolidNodeList1d_iterator", allow_subclassing=True)
        self.vector_of_SolidNodeList2d_iterator = addObject(mod, "vector_of_SolidNodeList2d_iterator", allow_subclassing=True)
        self.vector_of_SolidNodeList3d_iterator = addObject(mod, "vector_of_SolidNodeList3d_iterator", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateSolidFieldNamesBindings(self.SolidFieldNames)

        self.generateSolidNodeListBindings(self.SolidNodeList1d, 1)
        self.generateSolidNodeListBindings(self.SolidNodeList2d, 2)
        self.generateSolidNodeListBindings(self.SolidNodeList3d, 3)

        generateStdVectorBindings(self.vector_of_SolidNodeList1d, "Spheral::SolidMaterial::SolidNodeList1d*", "vector_of_SolidNodeList1d")
        generateStdVectorBindings(self.vector_of_SolidNodeList2d, "Spheral::SolidMaterial::SolidNodeList2d*", "vector_of_SolidNodeList2d")
        generateStdVectorBindings(self.vector_of_SolidNodeList3d, "Spheral::SolidMaterial::SolidNodeList3d*", "vector_of_SolidNodeList3d")
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

