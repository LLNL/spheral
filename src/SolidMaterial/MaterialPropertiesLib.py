#-------------------------------------------------------------------------------
# A python dictionary of material properties for various materials.
#-------------------------------------------------------------------------------
SpheralMaterialPropertiesLib = {

#-------------------------------------------------------------------------------
    "steel" : {
        "rho0" : 7.9,             # gm/cm^3
        "atomicWeight" : 55.350,  # gm/mol
        "Gruneisen" : {
            "C0"     : 0.457*1e6, # cm/sec
            "S1"     : 1.49,      # dimensionless
            "S2"     : 0.0,       # dimensionless
            "S3"     : 0.0,       # dimensionless
            "gamma0" : 1.93,      # dimensionless
            "b"      : 0.5,       # dimensionless (beta)
        },
        "kWeibull" : 6.52e5,      # cm^-3
        "mWeibull" : 2.63,        # dimensionless
    },

#-------------------------------------------------------------------------------
    "pumice" : {
        "rho0" : 2.327,         # gm/cm^3
        "atomicWeight" : 24.82, # gm/mol
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
        "atomicWeight" : 226.32,# gm/mol
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
        "atomicWeight" : 60.08, # gm/mol
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
    "gabbroic anorthosite 1pp melosh89" : {
        "rho0" : 2.940,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 71.0e10,       # (dyne/cm^2)
            "B"     : 75.0e10,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 487.0e10,     # erg/gm
            "epsLiquid" : 4.72e10,  # erg/gm
            "epsVapor"  : 18.2e10,  # erg/gm
        },
    },

#-------------------------------------------------------------------------------
    "gabbroic anorthosite hpp melosh89" : {
        "rho0" : 3.970,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.3,           # dimensionless
            "A"    : 240.0e10,      # (dyne/cm^2)
            "B"     : 130.0e10,     # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 1800.0e10,    # erg/gm
            "epsLiquid" : 3.19e10,  # erg/gm
            "epsVapor"  : 16.8e10,  # erg/gm
        },
    },

# 1 ----------------------------------------------------------------------------
    "granite" : {
        "rho0" : 2.680,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.50e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 1.00e27,       # cm^-3
        "mWeibull" : 6.2,           # dimensionless
    },

# 2 ----------------------------------------------------------------------------
    "basalt" : {
        "rho0" : 2.900,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.27e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 5.00e24,       # cm^-3
        "mWeibull" : 9.0,           # dimensionless
        "IvanovDamageModel" : {
            "epsfb" : 1e-4,         # dimensionless
            "B"     : 1e-12,        # 1/P (cm^2/dyne)
            "Pc"    : 3e9,          # dyne/cm^2
            "Yt"    : 10e7,         # dyne/cm^2, tensile failure stress
        },
        # "JohnsonCookStrength" : {
        #     "A"         : 244e7,    # dyne/cm^2
        #     "B"         : 488e7,    # dyne/cm^2
        #     "C"         : 0.02,     # dimensionless
        #     "C4"        : 0.0,      # dimensionless
        #     "m"         : 1.7,      # dimensionless
        #     "nhard"     : 0.5,      # dimensionless
        #     "epsdot0"   : 1.0,
        #     "epsdotmin" : 0.0,
        #     "Tmelt"     : 0.0,
        #     "Troom"     : 0.0,
        #     "mu0"       : 0.0,
        #     "shearModulusScaling" : False,
        # },
    },

# 3 ----------------------------------------------------------------------------
    "aluminum" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 24.032,# gm/mol
        "Gruneisen" : {
            "C0"     : 0.5386*1e6,# cm/sec
            "S1"     : 1.339,     # dimensionless
            "S2"     : 0.0,       # dimensionless
            "S3"     : 0.0,       # dimensionless
            "gamma0" : 1.97,      # dimensionless
            "b"      : 0.48,      # dimensionless (beta)
        },
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
        "mu0": 2.65e11,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2

        # Weibull parameters lifted from plots in DOI: 10.1007/s11661-010-0364-6
        "kWeibull" : 7e11/2.63**2.55, # cm^-3
        "mWeibull" : 2.55,            # dimensionless
    },

# 3 ----------------------------------------------------------------------------
    "aluminum melosh89" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 24.032,# gm/mol
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.63,          # dimensionless
            "A"    : 7.52e11,       # (dyne/cm^2)
            "B"     : 6.50e11,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 5.00e10,      # erg/gm
            "epsLiquid" : 3.00e10,  # erg/gm
            "epsVapor"  : 1.39e11,  # erg/gm
        },
        "mu0": 2.65e11,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2

        # Weibull parameters lifted from plots in DOI: 10.1007/s11661-010-0364-6
        "kWeibull" : 7e11/2.63**2.55, # cm^-3
        "mWeibull" : 2.55,            # dimensionless
    },

# 4 ----------------------------------------------------------------------------
    "copper" : {
        "rho0" : 8.900,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 0.0,                 # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
    },

# 5 ----------------------------------------------------------------------------
    "iron 130pt" : {
        "rho0" : 7.860,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 0.0,                 # dyne/cm^2
        "Y0" : 6.00e9,              # dyne/cm^2
    },

# 5 ----------------------------------------------------------------------------
    "iron melosh89" : {
        "rho0" : 7.800,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
        "Tillotson" : {
            "a"    : 0.5,           # dimensionless
            "b"    : 1.5,           # dimensionless
            "A"    : 1.28e12,       # (dyne/cm^2)
            "B"     : 1.05e12,      # (dyne/cm^2)
            "alpha" : 5.0,          # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 9.50e10,      # erg/gm
            "epsLiquid" : 2.40e10,  # erg/gm
            "epsVapor"  : 8.67e10,  # erg/gm
        },
        "mu0": 0.0,                 # dyne/cm^2
        "Y0" : 6.00e9,              # dyne/cm^2
    },

# 6 ----------------------------------------------------------------------------
    "lucite" : {
        "rho0" : 1.180,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 7.30e8,              # dyne/cm^2
        "Y0" : 1.00e9,              # dyne/cm^2
    },

# 7 ----------------------------------------------------------------------------
    "limestone" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.50e11,             # dyne/cm^2
        "Y0" : 1.00e33,             # dyne/cm^2
        "kWeibull" : 6.00e42,       # cm^-3
        "mWeibull" : 12.8,          # dimensionless
    },

# 8 ----------------------------------------------------------------------------
    "halite" : {
        "rho0" : 2.160,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 3.00e11,             # dyne/cm^2
        "Y0" : 9.00e8,              # dyne/cm^2
        "kWeibull" : 3.00e38,       # cm^-3
        "mWeibull" : 8.7,           # dimensionless
    },

# 9 ----------------------------------------------------------------------------
    "oil shale" : {
        "rho0" : 2.300,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 1.40e11,             # dyne/cm^2
        "Y0" : 1.00e33,             # dyne/cm^2
        "kWeibull" : 2.00e21,       # cm^-3
        "mWeibull" : 8.1,           # dimensionless
    },

#10 ----------------------------------------------------------------------------
    "wet tuff" : {
        "rho0" : 1.970,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 4.00e10,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
    },

#11 ----------------------------------------------------------------------------
    "dry tuff" : {
        "rho0" : 1.700,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.00e11,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
    },

#12 ----------------------------------------------------------------------------
    "alluvium" : {
        "rho0" : 2.600,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 0.0,                 # dyne/cm^2
        "Y0" : 1.00e9,              # dyne/cm^2
    },

#13 ----------------------------------------------------------------------------
    "anorthosite 1pp" : {
        "rho0" : 2.867,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 8.30e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 1.39e12,       # cm^-3
        "mWeibull" : 3.0,           # dimensionless
    },

#14 ----------------------------------------------------------------------------
    "anorthosite hpp" : {
        "rho0" : 3.970,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 8.30e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 5.00e22,       # cm^-3
        "mWeibull" : 9.1,           # dimensionless
    },

#15 ----------------------------------------------------------------------------
    "andesite" : {
        "rho0" : 2.700,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 4.00e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 5.00e22,       # cm^-3
        "mWeibull" : 8.5,           # dimensionless
    },

#16 ----------------------------------------------------------------------------
    "water" : {
        "rho0" : 0.998,         # gm/cm^3
        "atomicWeight" : 18.015,# gm/mol
        "Tillotson" : {
            "a"    : 0.7,           # dimensionless
            "b"    : 0.15,          # dimensionless
            "A"    : 2.18e10,       # (dyne/cm^2)
            "B"     : 13.25e10,     # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 7.00e10,      # erg/gm
            "epsLiquid" : 4.19e9,   # erg/gm
            "epsVapor"  : 2.69e10,  # erg/gm
        },
        "mu0": 0.0,                 # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
    },

#17 ----------------------------------------------------------------------------
    "pure ice" : {
        "rho0" : 0.917,         # gm/cm^3
        "atomicWeight" : 18.015,# gm/mol
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
        "mu0": 2.80e10,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
        "kWeibull" : 1.42e32,       # cm^-3
        "mWeibull" : 9.59,           # dimensionless
    },

#18 ----------------------------------------------------------------------------
    "5% silicate ice" : {
        "rho0" : 0.948,         # gm/cm^3
        "atomicWeight" : 18.015,# gm/mol
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
        "mu0": 2.80e10,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
        "kWeibull" : 5.60e37,       # cm^-3
        "mWeibull" : 9.4,           # dimensionless
    },

#19 ----------------------------------------------------------------------------
    "30% silicate ice"  : {
        "rho0" : 1.141,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.80e10,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
        "kWeibull" : 5.60e38,       # cm^-3
        "mWeibull" : 9.4,           # dimensionless
    },

#19 ----------------------------------------------------------------------------
    "ice melosh89"  : {
        "rho0" : 1.293,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
        "Tillotson" : {
            "a"    : 0.3,           # dimensionless
            "b"    : 0.1,           # dimensionless
            "A"    : 10.7e10,       # (dyne/cm^2)
            "B"     : 65.0e10,      # (dyne/cm^2)
            "alpha" : 10.0,         # dimensionless
            "beta"  : 5.0,          # dimensionless
            "eps0"  : 10.0e10,      # erg/gm
            "epsLiquid" : 7.73e9,   # erg/gm
            "epsVapor"  : 3.04e10,  # erg/gm
        },
        "mu0": 2.80e10,             # dyne/cm^2
        "Y0" : 1.00e10,             # dyne/cm^2
        "kWeibull" : 5.60e38,       # cm^-3
        "mWeibull" : 9.4,           # dimensionless
    },

#20 ----------------------------------------------------------------------------
    "special" : {
        "rho0" : 1.130,         # gm/cm^3
        "atomicWeight" : 60.08, # gm/mol
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
        "mu0": 2.27e11,             # dyne/cm^2
        "Y0" : 3.50e10,             # dyne/cm^2
        "kWeibull" : 5.00e24,       # cm^-3
        "mWeibull" : 9.0,           # dimensionless
    },

}
