#-------------------------------------------------------------------------------
# makeTillotsonEquationOfState
#
# Provides canned values for a variety of materials to build Tillotson EOS's.
#-------------------------------------------------------------------------------
from SpheralModules.Spheral.SolidMaterial import *
from SpheralModules.Spheral.Material import PhysicalConstants
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

#-------------------------------------------------------------------------------
# The base units for parameters in this file.
#-------------------------------------------------------------------------------
CGS = PhysicalConstants(0.01,    # Length in m
                        0.001,   # Mass in kg
                        1.0)     # Time in sec

#-------------------------------------------------------------------------------
# The generic factory function, where you pass in the dimension specific 
# Tillotson constructor.
# This one is for internal use only -- people will actually call the dimension
# specific front-ends at the end of this script.
#-------------------------------------------------------------------------------
def _TillotsonFactory(materialName,        # Label for the material in _TillotsonParams
                      etamin,              # Lower bound for rho/rho0
                      etamax,              # Upper bound for rho/rho0
                      units,               # Units the user wants to work in
                      externalPressure,    # Optional external pressure
                      minimumPressure,     # Optional minimum pressure
                      maximumPressure,     # Optional maximum pressure
                      TillConstructor):    # The actual dimension specific Tillotson constructor

    # Check that the caller specified a valid material label.
    mat = materialName.lower()
    if mat not in SpheralMaterialPropertiesLib:
        raise ValueError, "You must specify one of %s" % str(SpheralMaterialPropertiesLib.keys())
    if "Tillotson" not in SpheralMaterialPropertiesLib[mat]:
        raise ValueError, "The material %s does not provide Tillotson paramters." % materialName

    # Extract the parameters for this material.
    params = dict(SpheralMaterialPropertiesLib[mat]["Tillotson"])
    
    # Figure out the conversions to the requested units.
    lconv = CGS.unitLengthMeters / units.unitLengthMeters
    mconv = CGS.unitMassKg / units.unitMassKg
    tconv = CGS.unitTimeSec / units.unitTimeSec
    rhoConv = mconv/(lconv*lconv*lconv)
    Pconv = mconv/(lconv*tconv*tconv)
    specificEconv = (lconv/tconv)**2

    # Build the arguments for constructing the Tillotson.
    args = [SpheralMaterialPropertiesLib[mat]["rho0"] * rhoConv,
            etamin,
            etamax,
            params["a"],
            params["b"],
            params["A"] * Pconv,
            params["B"] * Pconv,
            params["alpha"],
            params["beta"],
            params["eps0"] * specificEconv,
            params["epsLiquid"] * specificEconv,
            params["epsVapor"] * specificEconv,
            SpheralMaterialPropertiesLib[mat]["atomicWeight"],
            units]
    if externalPressure:
        args.append(externalPressure)
    if minimumPressure:
        args.append(minimumPressure)
    if maximumPressure:
        args.append(maximumPressure)
    
    # Return the EOS.
    return TillConstructor(*tuple(args))

#-------------------------------------------------------------------------------
# Create the dimension specific Tillotson factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
# 1D
def makeTillotsonEquationOfState1d(materialName,
                                   etamin,
                                   etamax,
                                   units,
                                   externalPressure = None,
                                   minimumPressure = None,
                                   maximumPressure = None):
    return _TillotsonFactory(materialName,
                             etamin,
                             etamax,
                             units,
                             externalPressure,
                             minimumPressure,
                             maximumPressure,
                             TillotsonEquationOfState1d)

# 2D
def makeTillotsonEquationOfState2d(materialName,
                                   etamin,
                                   etamax,
                                   units,
                                   externalPressure = None,
                                   minimumPressure = None,
                                   maximumPressure = None):
    return _TillotsonFactory(materialName,
                             etamin,
                             etamax,
                             units,
                             externalPressure,
                             minimumPressure,
                             maximumPressure,
                             TillotsonEquationOfState2d)

# 3D
def makeTillotsonEquationOfState3d(materialName,
                                   etamin,
                                   etamax,
                                   units,
                                   externalPressure = None,
                                   minimumPressure = None,
                                   maximumPressure = None):
    return _TillotsonFactory(materialName,
                             etamin,
                             etamax,
                             units,
                             externalPressure,
                             minimumPressure,
                             maximumPressure,
                             TillotsonEquationOfState3d)
