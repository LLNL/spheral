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
        space = mod.add_cpp_namespace("Spheral")

        self.SolidFieldNames = addObject(space, "SolidFieldNames")

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateSolidFieldNamesBindings(self.SolidFieldNames)

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidFieldNames).
    #---------------------------------------------------------------------------
    def generateSolidFieldNamesBindings(self, x):
        x.add_static_attribute("deviatoricStress", "std::string",  is_const=True)
        x.add_static_attribute("deviatoricStressTT", "std::string",  is_const=True)
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
        x.add_static_attribute("flaws", "std::string",  is_const=True)
        x.add_static_attribute("effectiveFlaws", "std::string",  is_const=True)
        x.add_static_attribute("porosityAlpha", "std::string",  is_const=True)
        x.add_static_attribute("porosityStrain", "std::string",  is_const=True)
        x.add_static_attribute("fragmentIDs", "std::string",  is_const=True)
        x.add_static_attribute("particleTypes", "std::string",  is_const=True)
        x.add_static_attribute("meltSpecificEnergy", "std::string",  is_const=True)
        return
