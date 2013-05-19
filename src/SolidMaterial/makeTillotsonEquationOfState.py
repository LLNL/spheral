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
_TillotsonParams = {"pumice"           : {"rho0" : 2.327,         # gm/cm^3
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
                    "nylon"            : {"rho0" : 1.185,         # gm/cm^3
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
                    "glass"            : {"rho0" : 2.560,         # gm/cm^3
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
                    "granite"          : {"rho0" : 2.680,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.3,           # dimensionless
                                          "A"    : 1.80e11,       # (dyne/cm^2)
                                          "B"     : 1.80e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.60e11,      # erg/gm
                                          "epsLiquid" : 3.50e10,  # erg/gm
                                          "epsVapor"  : 1.80e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "basalt"           : {"rho0" : 2.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.5,           # dimensionless
                                          "A"    : 2.67e11,       # (dyne/cm^2)
                                          "B"     : 2.67e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 4.87e12,      # erg/gm
                                          "epsLiquid" : 4.72e10,  # erg/gm
                                          "epsVapor"  : 1.82e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "aluminum"         : {"rho0" : 2.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.63,          # dimensionless
                                          "A"    : 7.52e11,       # (dyne/cm^2)
                                          "B"     : 6.50e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 5.00e10,      # erg/gm
                                          "epsLiquid" : 2.70e10,  # erg/gm
                                          "epsVapor"  : 1.41e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "copper"           : {"rho0" : 8.900,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.5,           # dimensionless
                                          "A"    : 1.39e12,       # (dyne/cm^2)
                                          "B"     : 1.10e12,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 3.25e11,      # erg/gm
                                          "epsLiquid" : 1.35e10,  # erg/gm
                                          "epsVapor"  : 3.00e10,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "iron 130pt"       : {"rho0" : 7.860,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.5,           # dimensionless
                                          "A"    : 1.28e12,       # (dyne/cm^2)
                                          "B"     : 1.05e12,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 9.50e10,      # erg/gm
                                          "epsLiquid" : 1.42e10,  # erg/gm
                                          "epsVapor"  : 8.45e10,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "lucite"           : {"rho0" : 1.180,         # gm/cm^3
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
                    "limestone"        : {"rho0" : 2.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 0.6,           # dimensionless
                                          "A"    : 4.00e11,       # (dyne/cm^2)
                                          "B"     : 6.70e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.00e11,      # erg/gm
                                          "epsLiquid" : 2.50e10,  # erg/gm
                                          "epsVapor"  : 1.40e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "halite"           : {"rho0" : 2.160,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 0.6,           # dimensionless
                                          "A"    : 2.50e11,       # (dyne/cm^2)
                                          "B"     : 3.00e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 5.00e10,      # erg/gm
                                          "epsLiquid" : 2.00e10,  # erg/gm
                                          "epsVapor"  : 1.50e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "oil shale"        : {"rho0" : 2.300,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.0,           # dimensionless
                                          "A"    : 2.80e11,       # (dyne/cm^2)
                                          "B"     : 1.10e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.10e11,      # erg/gm
                                          "epsLiquid" : 3.20e10,  # erg/gm
                                          "epsVapor"  : 1.60e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "wet tuff"         : {"rho0" : 1.970,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.3,           # dimensionless
                                          "A"    : 1.00e11,       # (dyne/cm^2)
                                          "B"     : 6.00e10,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.10e11,      # erg/gm
                                          "epsLiquid" : 3.20e10,  # erg/gm
                                          "epsVapor"  : 1.60e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "dry tuff"         : {"rho0" : 1.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.3,           # dimensionless
                                          "A"    : 4.50e10,       # (dyne/cm^2)
                                          "B"     : 3.00e10,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 6.00e10,      # erg/gm
                                          "epsLiquid" : 3.50e10,  # erg/gm
                                          "epsVapor"  : 1.80e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "alluvium"         : {"rho0" : 2.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 0.8,           # dimensionless
                                          "A"    : 3.00e11,       # (dyne/cm^2)
                                          "B"     : 1.00e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 6.00e10,      # erg/gm
                                          "epsLiquid" : 3.50e10,  # erg/gm
                                          "epsVapor"  : 1.80e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "anorthosite 1pp"  : {"rho0" : 2.867,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.5,           # dimensionless
                                          "A"    : 7.10e11,       # (dyne/cm^2)
                                          "B"     : 7.50e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 4.87e12,      # erg/gm
                                          "epsLiquid" : 4.72e10,  # erg/gm
                                          "epsVapor"  : 1.82e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "anorthosite hpp"  : {"rho0" : 3.970,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.3,           # dimensionless
                                          "A"    : 2.40e12,       # (dyne/cm^2)
                                          "B"     : 1.30e12,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.80e13,      # erg/gm
                                          "epsLiquid" : 3.19e10,  # erg/gm
                                          "epsVapor"  : 1.68e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "andesite"         : {"rho0" : 2.700,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.3,           # dimensionless
                                          "A"    : 1.80e11,       # (dyne/cm^2)
                                          "B"     : 1.80e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.60e11,      # erg/gm
                                          "epsLiquid" : 3.50e10,  # erg/gm
                                          "epsVapor"  : 1.80e11,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "water"            : {"rho0" : 0.998,         # gm/cm^3
                                          "a"    : 0.7,           # dimensionless
                                          "b"    : 0.15,          # dimensionless
                                          "A"    : 2.18e10,       # (dyne/cm^2)
                                          "B"     : 1.33e11,      # (dyne/cm^2)
                                          "alpha" : 10.0,         # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 7.00e10,      # erg/gm
                                          "epsLiquid" : 4.19e9,   # erg/gm
                                          "epsVapor"  : 2.69e10,  # erg/gm
                                          "atomicWeight" : 18.015,# dimensionless
                                          },
                    "pure ice"         : {"rho0" : 0.917,         # gm/cm^3
                                          "a"    : 0.3,           # dimensionless
                                          "b"    : 0.1,           # dimensionless
                                          "A"    : 9.47e10,       # (dyne/cm^2)
                                          "B"     : 1.33e11,      # (dyne/cm^2)
                                          "alpha" : 10.0,         # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.00e11,      # erg/gm
                                          "epsLiquid" : 7.73e9,   # erg/gm
                                          "epsVapor"  : 3.04e10,  # erg/gm
                                          "atomicWeight" : 18.015,# dimensionless
                                          },
                    "5% silicate ice"  : {"rho0" : 0.948,         # gm/cm^3
                                          "a"    : 0.3,           # dimensionless
                                          "b"    : 0.1,           # dimensionless
                                          "A"    : 6.50e10,       # (dyne/cm^2)
                                          "B"     : 1.33e11,      # (dyne/cm^2)
                                          "alpha" : 10.0,         # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 1.00e11,      # erg/gm
                                          "epsLiquid" : 7.73e9,   # erg/gm
                                          "epsVapor"  : 3.04e10,  # erg/gm
                                          "atomicWeight" : 60.08, # dimensionless
                                          },
                    "30% silicate ice"  : {"rho0" : 1.141,         # gm/cm^3
                                           "a"    : 0.3,           # dimensionless
                                           "b"    : 0.1,           # dimensionless
                                           "A"    : 8.44e10,       # (dyne/cm^2)
                                           "B"     : 1.33e11,      # (dyne/cm^2)
                                           "alpha" : 10.0,         # dimensionless
                                           "beta"  : 5.0,          # dimensionless
                                           "eps0"  : 1.00e11,      # erg/gm
                                           "epsLiquid" : 7.73e9,   # erg/gm
                                           "epsVapor"  : 3.04e10,  # erg/gm
                                           "atomicWeight" : 60.08, # dimensionless
                                           },
                    "special"          : {"rho0" : 1.130,         # gm/cm^3
                                          "a"    : 0.5,           # dimensionless
                                          "b"    : 1.5,           # dimensionless
                                          "A"    : 2.67e11,       # (dyne/cm^2)
                                          "B"     : 2.67e11,      # (dyne/cm^2)
                                          "alpha" : 5.0,          # dimensionless
                                          "beta"  : 5.0,          # dimensionless
                                          "eps0"  : 4.87e12,      # erg/gm
                                          "epsLiquid" : 4.72e10,  # erg/gm
                                          "epsVapor"  : 1.82e11,  # erg/gm
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
    mat = materialName.lower()
    if mat not in _TillotsonParams:
        raise ValueError, "You must specify one of %s" % str(_TillotsonParams.keys())

    # Extract the parameters for this material.
    params = dict(_TillotsonParams[mat])
    
    # Figure out the conversions to the requested units.
    lconv = CGS.unitLengthMeters / units.unitLengthMeters
    mconv = CGS.unitMassKg / units.unitMassKg
    tconv = CGS.unitTimeSec / units.unitTimeSec
    rhoConv = mconv/(lconv*lconv*lconv)
    Pconv = mconv/(lconv*tconv*tconv)
    specificEconv = (lconv/tconv)**2

    # Build the arguments for constructing the Tillotson.
    args = [params["rho0"] * rhoConv,
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
