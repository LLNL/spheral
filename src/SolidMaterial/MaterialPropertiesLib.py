#-------------------------------------------------------------------------------
# A python dictionary of material properties for various materials.
#-------------------------------------------------------------------------------
SpheralMaterialPropertiesLib = {

#-------------------------------------------------------------------------------
    "pumice" : {
        "rho0" : 2.327,         # gm/cm^3
        "atomicWeight" : 24.82, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 2.67e11,       # (dyne/cm^2)
            "B"     : 2.67e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 4.87e12,      # erg/gm
            "epsLiquid" : 4.72e10,  # erg/gm
            "epsVapor"  : 1.82e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "nylon" : {
        "rho0" : 1.185,         # gm/cm^3
        "atomicWeight" : 226.32,# dimensionless
        "Tillotson" : {
            "a"    : 0.6,           # dimensionless
            "b"    : 2.0,           # dimensionless
            "A"    : 1.01e11,       # (dyne/cm^2)
            "B"     : 3.38e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 7.00e10,      # erg/gm
            "epsLiquid" : 2.00e10,  # erg/gm
            "epsVapor"  : 2.40e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "glass" : {
        "rho0" : 2.560,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.6,           # dimensionless
            "b"    : 2.0,           # dimensionless
            "A"    : 1.01e11,       # (dyne/cm^2)
            "B"     : 3.38e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 7.00e10,      # erg/gm
            "epsLiquid" : 2.00e10,  # erg/gm
            "epsVapor"  : 2.40e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "granite" : {
        "rho0" : 2.680,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 1.80e11,       # (dyne/cm^2)
            "B"     : 1.80e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.60e11,      # erg/gm
            "epsLiquid" : 3.50e10,  # erg/gm
            "epsVapor"  : 1.80e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "basalt" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 2.67e11,       # (dyne/cm^2)
            "B"     : 2.67e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 4.87e12,      # erg/gm
            "epsLiquid" : 4.72e10,  # erg/gm
            "epsVapor"  : 1.82e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "aluminum" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.63,          # dimensionless
            "A"    : 7.52e11,       # (dyne/cm^2)
            "B"     : 6.50e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 5.00e10,      # erg/gm
            "epsLiquid" : 2.70e10,  # erg/gm
            "epsVapor"  : 1.41e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "copper" : {
        "rho0" : 8.900,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 1.39e12,       # (dyne/cm^2)
            "B"     : 1.10e12,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 3.25e11,      # erg/gm
            "epsLiquid" : 1.35e10,  # erg/gm
            "epsVapor"  : 3.00e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "iron 130pt" : {
        "rho0" : 7.860,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 1.28e12,       # (dyne/cm^2)
            "B"     : 1.05e12,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 9.50e10,      # erg/gm
            "epsLiquid" : 1.42e10,  # erg/gm
            "epsVapor"  : 8.45e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "lucite" : {
        "rho0" : 1.180,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.6,           # dimensionless
            "b"    : 2.0,           # dimensionless
            "A"    : 1.01e11,       # (dyne/cm^2)
            "B"     : 3.38e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 7.00e10,      # erg/gm
            "epsLiquid" : 2.00e10,  # erg/gm
            "epsVapor"  : 2.40e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "limestone" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 0.6,           # dimensionless
            "A"    : 4.00e11,       # (dyne/cm^2)
            "B"     : 6.70e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.00e11,      # erg/gm
            "epsLiquid" : 2.50e10,  # erg/gm
            "epsVapor"  : 1.40e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "halite" : {
        "rho0" : 2.160,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 0.6,           # dimensionless
            "A"    : 2.50e11,       # (dyne/cm^2)
            "B"     : 3.00e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 5.00e10,      # erg/gm
            "epsLiquid" : 2.00e10,  # erg/gm
            "epsVapor"  : 1.50e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "oil shale" : {
        "rho0" : 2.300,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.0,           # dimensionless
            "A"    : 2.80e11,       # (dyne/cm^2)
            "B"     : 1.10e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.10e11,      # erg/gm
            "epsLiquid" : 3.20e10,  # erg/gm
            "epsVapor"  : 1.60e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "wet tuff" : {
        "rho0" : 1.970,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 1.00e11,       # (dyne/cm^2)
            "B"     : 6.00e10,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.10e11,      # erg/gm
            "epsLiquid" : 3.20e10,  # erg/gm
            "epsVapor"  : 1.60e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "dry tuff" : {
        "rho0" : 1.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 4.50e10,       # (dyne/cm^2)
            "B"     : 3.00e10,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 6.00e10,      # erg/gm
            "epsLiquid" : 3.50e10,  # erg/gm
            "epsVapor"  : 1.80e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "alluvium" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 0.8,           # dimensionless
            "A"    : 3.00e11,       # (dyne/cm^2)
            "B"     : 1.00e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 6.00e10,      # erg/gm
            "epsLiquid" : 3.50e10,  # erg/gm
            "epsVapor"  : 1.80e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "anorthosite 1pp" : {
        "rho0" : 2.867,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 7.10e11,       # (dyne/cm^2)
            "B"     : 7.50e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 4.87e12,      # erg/gm
            "epsLiquid" : 4.72e10,  # erg/gm
            "epsVapor"  : 1.82e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "anorthosite hpp" : {
        "rho0" : 3.970,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 2.40e12,       # (dyne/cm^2)
            "B"     : 1.30e12,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.80e13,      # erg/gm
            "epsLiquid" : 3.19e10,  # erg/gm
            "epsVapor"  : 1.68e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "andesite" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 1.80e11,       # (dyne/cm^2)
            "B"     : 1.80e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.60e11,      # erg/gm
            "epsLiquid" : 3.50e10,  # erg/gm
            "epsVapor"  : 1.80e11,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "water" : {
        "rho0" : 0.998,         # gm/cm^3
        "atomicWeight" : 18.015,# dimensionless
        "Tillotson" : {
            "a"    : 0.7,           # dimensionless
            "b"    : 0.15,          # dimensionless
            "A"    : 2.18e10,       # (dyne/cm^2)
            "B"     : 1.33e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 7.00e10,      # erg/gm
            "epsLiquid" : 4.19e9,   # erg/gm
            "epsVapor"  : 2.69e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "pure ice" : {
        "rho0" : 0.917,         # gm/cm^3
        "atomicWeight" : 18.015,# dimensionless
        "Tillotson" : {
            "a"    : 0.3,           # dimensionless
            "b"    : 0.1,           # dimensionless
            "A"    : 9.47e10,       # (dyne/cm^2)
            "B"     : 1.33e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.00e11,      # erg/gm
            "epsLiquid" : 7.73e9,   # erg/gm
            "epsVapor"  : 3.04e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "5% silicate ice" : {
        "rho0" : 0.948,         # gm/cm^3
        "atomicWeight" : 18.015,# dimensionless
        "Tillotson" : {
            "a"    : 0.3,           # dimensionless
            "b"    : 0.1,           # dimensionless
            "A"    : 6.50e10,       # (dyne/cm^2)
            "B"     : 1.33e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.00e11,      # erg/gm
            "epsLiquid" : 7.73e9,   # erg/gm
            "epsVapor"  : 3.04e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "30% silicate ice"  : {
        "rho0" : 1.141,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.3,           # dimensionless
            "b"    : 0.1,           # dimensionless
            "A"    : 8.44e10,       # (dyne/cm^2)
            "B"     : 1.33e11,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1.00e11,      # erg/gm
            "epsLiquid" : 7.73e9,   # erg/gm
            "epsVapor"  : 3.04e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "special" : {
        "rho0" : 1.130,         # gm/cm^3
        "atomicWeight" : 60.08, # dimensionless
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 2.67e11,       # (dyne/cm^2)
            "B"     : 2.67e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 4.87e12,      # erg/gm
            "epsLiquid" : 4.72e10,  # erg/gm
            "epsVapor"  : 1.82e11,  # erg/gm
        },
    },

}
