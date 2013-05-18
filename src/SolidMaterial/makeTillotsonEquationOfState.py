#-------------------------------------------------------------------------------
# makeTillotsonEquationOfState
#
# Provides canned values for a variety of materials to build Tillotson EOS's.
#-------------------------------------------------------------------------------
from SpheralModules.Spheral.SolidMaterial import *
from SpheralModules.Spheral.Material import PhysicalConstants

#-------------------------------------------------------------------------------
# A dictionary to provide the per material values.
# All values here expressed in CGS units.
#-------------------------------------------------------------------------------
_TillotsonParams = {"Pumice" : {"rho0" : 2.327,         # gm/cm^3
                                "a"    : 0.5,           # dimensionless
                                "b"    : 1.5,           # dimensionless
                                "A"    : 2.67e11,       # (dyne/cm^2)
                                "B"     : 2.67e11,      # (dyne/cm^2)
                                "alpha" : 5.0,          # dimensionless
                                "beta"  : 5.0,          # dimensionless
                                "eps0"  : 4.87e12,      # erg/gm
                                "epsLiquid" : 4.72e10,  # erg/gm
                                "epsVapor"  : 1.82e11,  # erg/gm
                                "atomicWeight" : 24.82, # dimensionless
                                },
                    "Nylon"  : {"rho0" : 1.185,         # gm/cm^3
                                "a"    : 0.6,           # dimensionless
                                "b"    : 2.0,           # dimensionless
                                "A"    : 1.01e11,       # (dyne/cm^2)
                                "B"     : 3.38e11,      # (dyne/cm^2)
                                "alpha" : 10.0,         # dimensionless
                                "beta"  : 5.0,          # dimensionless
                                "eps0"  : 7.00e10,      # erg/gm
                                "epsLiquid" : 2.00e10,  # erg/gm
                                "epsVapor"  : 2.40e10,  # erg/gm
                                "atomicWeight" : 226.32,# dimensionless
                                },
                    "Glass"  : {"rho0" : 2.560,         # gm/cm^3
                                "a"    : 0.6,           # dimensionless
                                "b"    : 2.0,           # dimensionless
                                "A"    : 1.01e11,       # (dyne/cm^2)
                                "B"     : 3.38e11,      # (dyne/cm^2)
                                "alpha" : 10.0,         # dimensionless
                                "beta"  : 5.0,          # dimensionless
                                "eps0"  : 7.00e10,      # erg/gm
                                "epsLiquid" : 2.00e10,  # erg/gm
                                "epsVapor"  : 2.40e10,  # erg/gm
                                "atomicWeight" : 60.08, # dimensionless
                                },
                    }

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
    if materialName not in _TillotsonParams:
        raise ValueError, "You must specify one of %s" % str(_TillotsonParams.keys())

    # Extract the parameters for this material.
    params = dict(_TillotsonParams[materialName])
    
    # Figure out the conversions to the requested units.
    lconv = units.unitLengthMeters / CGS.unitLengthMeters
    mconv = units.unitMassKg / CGS.unitMassKg
    tconv = units.unitTimeSec / CGS.unitTimeSec
    rhoConv = mconv/(lconv*lconv*lconv)
    Pconv = mconv/(lconv*tconv*tconv)
    specificEconv = (lconv/tconv)**2

    # Build the arguments for constructing the Tillotson.
    args = [params["rho0"] / rhoConv,
            etamin,
            etamax,
            params["a"],
            params["b"],
            params["A"] / Pconv,
            params["B"] / Pconv,
            params["alpha"],
            params["beta"],
            params["eps0"] / specificEconv,
            params["epsLiquid"] / specificEconv,
            params["epsVapor"] / specificEconv,
            params["atomicWeight"],
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
