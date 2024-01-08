#-------------------------------------------------------------------------------
# ShadowTillotsonEquationOfState
#
# Provides convenient constructors for the Tillotson using the canned values
# in MaterialPropertiesLib.py.
#-------------------------------------------------------------------------------
import types
from spheralDimensions import spheralDimensions
from CaptureStdout import helpString
from SpheralCompiledPackages import PhysicalConstants, PressureFloor, ZeroPressure
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

dims = spheralDimensions()
for dim in dims:
    exec("""
from SpheralCompiledPackages import TillotsonEquationOfState%(dim)sd as RealTillotsonEquationOfState%(dim)sd
""" % {"dim" : dim})

#-------------------------------------------------------------------------------
# Define a string providing the help for building a TillotsonEquationOfState.
#-------------------------------------------------------------------------------
expectedUsageString = """
TillotsonEquationOfState can be constructed one of two ways:

1.  Using canned material values stored in our internal data base.  Expected arguments:
        materialName          : Label for the material in data base                               
        etamin                : Lower bound for rho/rho0                                          
        etamax                : Upper bound for rho/rho0                                          
        units                 : Units the user wants to work in                                   
        etamin_solid          : Optional more restrictive lower bound for rho/rho0 in solid phase 
        etamax_solid          : Optional more restrictive upper bound for rho/rho0 in solid phase 
        externalPressure      : Optional external pressure                                        
        minimumPressure       : Optional minimum pressure                                         
        maximumPressure       : Optional maximum pressure                                         
        minimumPressureDamage : Optional minimum pressure in damage (default 0.0)
        minPressureType       : Optional behavior at minimumPressure (PressureFloor, ZeroPressure)

2.  You can directly set all the Tillotson parameters explicitly, as
        referenceDensity      : reference material mass density             
        etamin                : minimum allowed ratio rho/rho0              
        etamax                : maximum allowed ratio rho/rho0              
        etamin_solid          : minimum allowed ratio rho/rho0 (solid phase)
        etamax_solid          : maximum allowed ratio rho/rho0 (solid phase)
        a
        b
        A
        B
        alpha
        beta
        eps0
        epsLiquid
        epsVapor
        atomicWeight
        units
        externalPressure      : Optional external pressure                         
        minimumPressure       : Optional minimum pressure                          
        maximumPressure       : Optional maximum pressure                          
        minimumPressureDamage : Optional minimum pressure in damage (default 0.0)
        minPressureType       : Optional behavior at minimumPressure (PressureFloor, ZeroPressure)
"""

#-------------------------------------------------------------------------------
# The generic factory function, where you pass in the dimension specific 
# Tillotson constructor.
# This one is for internal use only -- people will actually call the dimension
# specific front-ends at the end of this script.
#-------------------------------------------------------------------------------
def _TillotsonFactory(*args, 
                      **kwargs):

    # The calling routine must provide the appropriate C++ constructor.
    TillConstructor = kwargs["TillConstructor"]

    # The arguments that need to be passed to this method.
    expectedArgs = ["materialName", "etamin", "etamax", "units"]
    optionalKwArgs = {"etamin_solid"          : 0.0,   
                      "etamax_solid"          : 1e200, 
                      "externalPressure"      : 0.0,   
                      "minimumPressure"       : -1e200,
                      "maximumPressure"       :  1e200,
                      "minimumPressureDamage" : 0.0,
                      "minPressureType"       : PressureFloor}

    # The base units for parameters in this file.
    CGS = PhysicalConstants(0.01,    # Length in m
                            0.001,   # Mass in kg
                            1.0)     # Time in sec
    # What sort of information did the user pass in?
    if ("materialName" in kwargs or 
        len(args) > 0 and type(args[0]) is str):

        # It looks like the user is trying to use one of the libarary canned values.
        # Evaluate the arguments to the method.
        if (len(args) > len(expectedArgs) or 
            (len(args) + len(kwargs) < len(expectedArgs))): # insist on formal mandatory arguments 
            raise ValueError(expectedUsageString)
        
        # Check for any invalid keywords
        for arg in kwargs: # deal with optional args
            if arg not in (expectedArgs + list(optionalKwArgs.keys()) + ["TillConstructor"]):
                raise ValueError(expectedUsageString)

        # Set the arguments dictionary
        dargs = {expectedArgs[i] : args[i] for i in range(len(args))}          # Mandatory args
        dargs.update(optionalKwArgs)
        dargs.update(kwargs)
        ARGS = types.SimpleNamespace(**dargs)

        # Check that the caller specified a valid material label.
        mat = ARGS.materialName.lower()
        if mat not in SpheralMaterialPropertiesLib:
            raise ValueError("You must specify one of %s" % str(list(SpheralMaterialPropertiesLib.keys())))
        if "Tillotson" not in SpheralMaterialPropertiesLib[mat]:
            raise ValueError("The material %s does not provide Tillotson paramters." % materialName)

        # Extract the parameters for this material.
        params = types.SimpleNamespace(**dict(SpheralMaterialPropertiesLib[mat]["Tillotson"]))
    
        # Figure out the conversions to the requested units.
        lconv = CGS.unitLengthMeters / ARGS.units.unitLengthMeters
        mconv = CGS.unitMassKg / ARGS.units.unitMassKg
        tconv = CGS.unitTimeSec / ARGS.units.unitTimeSec
        rhoConv = mconv/(lconv*lconv*lconv)
        Pconv = mconv/(lconv*tconv*tconv)
        specificEconv = (lconv/tconv)**2

        # Build the arguments for constructing the Tillotson.
        passargs = [SpheralMaterialPropertiesLib[mat]["rho0"] * rhoConv,
                    ARGS.etamin,
                    ARGS.etamax,
                    ARGS.etamin_solid,
                    ARGS.etamax_solid,
                    params.a,
                    params.b,
                    params.A * Pconv,
                    params.B * Pconv,
                    params.alpha,
                    params.beta,
                    params.eps0 * specificEconv,
                    params.epsLiquid * specificEconv,
                    params.epsVapor * specificEconv,
                    SpheralMaterialPropertiesLib[mat]["atomicWeight"] * mconv,
                    ARGS.units]
        passargs.extend([ARGS.externalPressure,
                         ARGS.minimumPressure,
                         ARGS.maximumPressure,
                         ARGS.minimumPressureDamage,
                         ARGS.minPressureType])
        passkwargs = {}

    else:

        # Just pass through the arguments.
        passargs = args
        passkwargs = kwargs
        del passkwargs["TillConstructor"]
    
    # Return the EOS.
    return TillConstructor(*tuple(passargs), **passkwargs)

#-------------------------------------------------------------------------------
# Create the dimension specific Tillotson factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
for dim in dims:
    exec("""
def TillotsonEquationOfState%(dim)sd(*args, **kwargs):
    expectedUsageString
    kwargs["TillConstructor"] = RealTillotsonEquationOfState%(dim)sd
    return _TillotsonFactory(*args, **kwargs)

TillotsonEquationOfState%(dim)sd.__doc__ = expectedUsageString + "\\n\\n" + helpString(RealTillotsonEquationOfState%(dim)sd)
""" % {"dim" : dim})
