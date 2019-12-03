#-------------------------------------------------------------------------------
# ShadowConstantStrength
#
# Provides convenient constructors for the ConstantStrength model using the canned
# values in MaterialPropertiesLib.py.
#-------------------------------------------------------------------------------
from SpheralCompiledPackages import PhysicalConstants
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

for dim in dims:
    exec("""
from SpheralCompiledPackages import ConstantStrength%(dim)sd as RealConstantStrength%(dim)sd
""" % {"dim" : dim})

#-------------------------------------------------------------------------------
# Define a string providing the help for building a ConstantStrength.
#-------------------------------------------------------------------------------
expectedUsageString = """
ConstantStrength can be constructed one of two ways:

1.  Using canned material values stored in our internal data base.  Expected arguments:
        materialName        : Label for the material in data base

2.  You can directly set the strength parameters explicitly, as
        mu0                 : shear modulus
        Y0                  : yield strength for plastic yielding
"""

#-------------------------------------------------------------------------------
# The generic factory function, where you pass in the dimension specific 
# ConstantStrength constructor.
# This one is for internal use only -- people will actually call the dimension
# specific front-ends at the end of this script.
#-------------------------------------------------------------------------------
def _ConstantStrengthFactory(*args, 
                              **kwargs):

    # The calling routine must provide the appropriate C++ constructor.
    CSConstructor = kwargs["CSConstructor"]

    # The arguments that need to be passed to this method.
    expectedArgs = ["materialName", "units"]
    optionalKwArgs = {"mu0" : None,
                      "Y0" : None}

    # The base units for parameters in this file.
    cgs = PhysicalConstants(0.01,    # Length in m
                            0.001,   # Mass in kg
                            1.0)     # Time in sec

    # What sort of information did the user pass in?
    if ("materialName" in kwargs or 
        len(args) > 0 and type(args[0]) is str):

        # It looks like the user is trying to use one of the libarary canned values.
        # Evaluate the arguments to the method.
        if len(args) > 0:
            if len(args) != len(expectedArgs):
                raise ValueError, expectedUsageString
            for i in xrange(len(expectedArgs)):
                exec("%s = args[i]" % expectedArgs[i])
            for arg in optionalKwArgs:
                exec("%s = optionalKwArgs['%s']" % (arg, arg))
        else:
            for arg in kwargs:
                if arg not in (expectedArgs + optionalKwArgs.keys() + ["CSConstructor"]):
                    raise ValueError, expectedUsageString
                exec("%s = kwargs['%s']" % (arg, arg))
            for arg in optionalKwArgs:
                if arg not in kwargs:
                    exec("%s = optionalKwArgs['%s']" % (arg, arg))

        # Check that the caller specified a valid material label.
        mat = materialName.lower()
        if mat not in SpheralMaterialPropertiesLib:
            raise ValueError, "You must specify one of %s" % str(SpheralMaterialPropertiesLib.keys())
        if ("mu0" not in SpheralMaterialPropertiesLib[mat] or
            "Y0"  not in SpheralMaterialPropertiesLib[mat]):
            raise ValueError, "The material %s does not provide strength paramters." % materialName

        # Extract the parameters for this material.
        if mu0 is None:
            mu0 = SpheralMaterialPropertiesLib[mat]["mu0"]
        if Y0 is None:
            Y0 = SpheralMaterialPropertiesLib[mat]["Y0"]
    
        # Figure out the conversions to the requested units.
        lconv = cgs.unitLengthMeters / units.unitLengthMeters
        mconv = cgs.unitMassKg / units.unitMassKg
        tconv = cgs.unitTimeSec / units.unitTimeSec
        rhoConv = mconv/(lconv*lconv*lconv)
        Pconv = mconv/(lconv*tconv*tconv)
        specificEconv = (lconv/tconv)**2

        # Build the arguments for constructing the ConstantStrength.
        passargs = [mu0 * Pconv,
                    Y0 * Pconv]
        passkwargs = {}

    else:

        # Just pass through the arguments.
        passargs = args
        passkwargs = kwargs
        del passkwargs["CSConstructor"]
    
    # Return the EOS.
    return CSConstructor(*tuple(passargs), **passkwargs)

#-------------------------------------------------------------------------------
# Create the dimension specific ConstantStrength factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
for dim in dims:
    exec("""
def ConstantStrength%(dim)sd(*args, **kwargs):
    expectedUsageString
    kwargs["CSConstructor"] = RealConstantStrength%(dim)sd
    return _ConstantStrengthFactory(*args, **kwargs)

ConstantStrength%(dim)sd.__doc__ = expectedUsageString
""" % {"dim" : dim})
