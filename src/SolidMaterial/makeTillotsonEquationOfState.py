#-------------------------------------------------------------------------------
# makeTillotsonEquationOfState
#
# Provides canned values for a variety of materials to build Tillotson EOS's.
#-------------------------------------------------------------------------------
from SpheralModules.Spheral.SolidMaterial import *
from SpheralModules.Spheral import PhysicalConstants

# A dictionary to provide the per material values.
# All values here expressed in CGS units.
_TillotsonParams = {"Pumice" : {"rho0" : 2.327,         # gm/cm^3
                                "a"    : 0.5,           # dimensionless
                                "b"    : 1.5,           # dimensionless
                                "A"    : 2.67e11        # (dyne/cm^2)
                                "B"     : 2.67e11       # (dyne/cm^2)
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
                                "A"    : 1.01e11        # (dyne/cm^2)
                                "B"     : 3.38e11       # (dyne/cm^2)
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
                                "A"    : 1.01e11        # (dyne/cm^2)
                                "B"     : 3.38e11       # (dyne/cm^2)
                                "alpha" : 10.0,         # dimensionless
                                "beta"  : 5.0,          # dimensionless
                                "eps0"  : 7.00e10,      # erg/gm
                                "epsLiquid" : 2.00e10,  # erg/gm
                                "epsVapor"  : 2.40e10,  # erg/gm
                                "atomicWeight" : 60.08, # dimensionless
                                },
                    }

# The base units for parameters in this file.
CGS = PhysicalConstants(0.01,    # Length in m
                        0.001,   # Mass in kg
                        1.0)     # Time in sec


# The generic factory function, where you pass in the dimension specific Tillotson constructor.
def _TillotsonFactory(materialName,
                      units,
                      externalPressure,
                      minimumPressure,
                      maximumPressure,
                      TillConstructor):
    if materialName not in TillotsonParams:
        raise ValueError, "You must specify one of %s" % str(TillotsonParams.keys())
    params = dict(_TillotsonParams)
    if externalPressure:
        params["externalPressure"] = externalPressure
    if minimumPressure:
        params["minimumPressure"] = minimumPressure
    if maximumPressure:
        params["maximumPressure"] = maximumPressure
    return TillConstructor(**params)

# Create the dimension specific Tillotson factories.  These are the ones
# you actually use.
def makeTillotsonEquationOfState
