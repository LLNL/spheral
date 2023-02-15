import sys
import mpi
import io, contextlib
from math import *
from SpheralCompiledPackages import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Define help strings for our constructors.
#-------------------------------------------------------------------------------
expectedUsageString = """
IvanoviSALEDamageModel is constructed with the following arguments (any 
default values listed in parens):

        materialName            : (optional) label for the material in data base to
                                  lookup the Weibull parameters
        units                   : (optional) scale kWeibull/mWeibull lookup from
                                  material name to given units
        nodeList                : (required) the FluidNodeList this model should be 
                                  applied to
        kernel                  : (required) the interpolation kernel to use
        minPlasticFailure       : min plastic strain for shearing plastic failure to 
                                  start
        plasticFailurePressureSlope  : slope for plastic strain shearing failure law
        plasticFailurePressureOffset : intercept for plastic strain shearing failure 
                                       law
        tensileFailureStress    : threshold stress for tensile failure to start
        crackGrowthMultiplier   : (0.4) crack growth rate in units of longitudinal
                                  sound speed
        damageCouplingAlgorithm : (DirectDamage) how should damaged points couple
        criticalDamageThreshold : (4.0) prevent any nodes where Trace(D_i) exceeds 
                                  criticalDamageThreshold from setting timestep
        mask                    : (1 on all points) a field of flags: a node with 
                                  zero implies no flaws (and therefore no damage) 
                                  on that point
"""

#-------------------------------------------------------------------------------
# IvanoviSALEDamageModel
#-------------------------------------------------------------------------------
IvanoviSALEDamageModelGenString = """
class IvanoviSALEDamageModel%(dim)s(CXXIvanoviSALEDamageModel%(dim)s):
    '''%(help)s'''

    def __init__(self, *args_in, **kwargs):

        args = list(args_in)

        # The material library values are in CGS, so build a units object for 
        # conversions.
        cgs = PhysicalConstants(0.01,   # unit length in m
                                0.001,  # unit mass in kg
                                1.0)    # unit time in sec

        # Arguments needed to build the damage model.
        damage_kwargs = {"nodeList"                     : None,
                         "kernel"                       : None,
                         "minPlasticFailure"            : None,
                         "plasticFailurePressureSlope"  : None,
                         "plasticFailurePressureOffset" : None,
                         "tensileFailureStress"         : None,
                         "crackGrowthMultiplier"        : 0.4,
                         "damageCouplingAlgorithm"      : DirectDamage,
                         "criticalDamageThreshold"      : 4.0,
                         "mask"                         : None}

        # Extra arguments for our convenient constructor.
        convenient_kwargs = {"materialName"             : None,
                             "units"                    : None}

        # Check the input arguments.
        validKeys = damage_kwargs.keys() + convenient_kwargs.keys()
        for argname in kwargs:
            if not argname in validKeys:
                raise ValueError("ERROR: argument %%s not a valid option.\\n" %% argname +
                                 expectedUsageString)

        # Did the user try any convenient constructor operations?
        if ((len(args) > 0 and type(args[0]) == str) or
            "materialName" in kwargs):
            if len(args) > 0 and type(args[0]) == str:
                materialName = args[0]
                del args[0]
            else:
                materialName = kwargs["materialName"]
                del kwargs["materialName"]
            if not materialName in SpheralMaterialPropertiesLib:
                raise ValueError(("ERROR: material %%s is not in the library of material values.\\n" %% materialName) +
                                 expectedUsageString)
            matprops = SpheralMaterialPropertiesLib[materialName]
            if not ("IvanovDamageModel" in matprops):
                raise ValueError(("ERROR : material %%s does not provide the required values for the Ivanov damage model.\\n" %% materialName) + 
                                 expectedUsageString)
            damage_kwargs["minPlasticFailure"] = matprops["IvanovDamageModel"]["epsfb"]
            damage_kwargs["plasticFailurePressureSlope"] = matprops["IvanovDamageModel"]["B"]
            damage_kwargs["plasticFailurePressureOffset"] = matprops["IvanovDamageModel"]["Pc"]
            damage_kwargs["tensileFailureStress"] = matprops["IvanovDamageModel"]["Yt"]

            # Any attempt to specify units?
            units = None
            if "units" in kwargs:
                units = kwargs["units"]
                del kwargs["units"]
            elif len(args) > 1 and isinstance(args[1], PhysicalConstants):
                units = args[1]
                del args[1]
            if units:
                  lconv = cgs.unitLengthMeters / units.unitLengthMeters
                  mconv = cgs.unitMassKg / units.unitMassKg
                  tconv = cgs.unitTimeSec / units.unitTimeSec
                  Pconv = mconv/(lconv*tconv*tconv)
                  damage_kwargs["plasticFailurePressureSlope"] /= Pconv
                  damage_kwargs["plasticFailurePressureOffset"] *= Pconv
                  damage_kwargs["tensileFailureStress"] *= Pconv

        # Process remaining user arguments.
        kwarg_order = ["nodeList",
                       "kernel",
                       "minPlasticFailure",
                       "plasticFailurePressureSlope",
                       "plasticFailurePressureOffset",
                       "tensileFailureStress",
                       "crackGrowthMultiplier",
                       "damageCouplingAlgorithm",
                       "damageInCompression",
                       "criticalDamageThreshold",
                       "mask"]
        for iarg, argval in enumerate(args):
            damage_kwargs[kward_order[iarg]] = argval

        # Process any keyword arguments.  Note we already removed any deprecated keywords.
        for argname in kwargs:
            if argname in damage_kwargs:
                damage_kwargs[argname] = kwargs[argname]
            else:
                raise ValueError(("ERROR : unknown kwarg %%s.\\n" %% argname) + expectedUsageString)

        # If no mask was provided, deafult to all points active
        if damage_kwargs["mask"] is None:
            damage_kwargs["mask"] = IntField%(dim)s("damage mask", damage_kwargs["nodeList"], 1)

        # Build the damage model.
        CXXIvanoviSALEDamageModel%(dim)s.__init__(self, **damage_kwargs)
        return

"""

#-------------------------------------------------------------------------------
# Make 'em
#-------------------------------------------------------------------------------
for dim in dims:
    exec("from SpheralCompiledPackages import IvanoviSALEDamageModel%id as CXXIvanoviSALEDamageModel%id" % (dim, dim))

    # Capture the full class help string
    with contextlib.redirect_stdout(io.StringIO()) as ss:
        exec("help(CXXIvanoviSALEDamageModel%id)" % dim)
    class_help = ss.getvalue()

    exec(IvanoviSALEDamageModelGenString % {"dim": "%id" % dim,
                                            "help": (expectedUsageString + "\n\n Class help:\n\n" + class_help)})
