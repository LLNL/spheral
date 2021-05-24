import sys
from io import BytesIO as StringIO
from SpheralCompiledPackages import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Define help strings for our constructors.
#-------------------------------------------------------------------------------
expectedUsageString = """
ProbabilisticDamageModel is constructed with the following arguments (any 
default values listed in parens):

        materialName            : (optional) label for the material in data base to
                                  lookup the Weibull parameters
        units                   : (optional) scale kWeibull/mWeibull lookup from
                                  material name to given units
        nodeList                : (required) the FluidNodeList this model should be 
                                  applied to
        kernel                  : (required) the interpolation kernel to use
        kWeibull                : the "k" Weibull constant -- can be looked up from 
                                  materialName or provided
        mWeibull                : the "m" Weibull constant -- can be looked up from 
                                  materialName or provided
        seed                    : (48927592) random number seed for flaw generation
        minFlawsPerNode         : (1) the minimum number of flaws to seed on a point
        crackGrowthMultiplier   : (0.4) crack growth rate in units of longitudinal
                                  sound speed
        volumeMultiplier        : (1.0) Multiplies per node volume, useful for 
                                  reduced dimensional (1D or 2D) problems.
        damageCouplingAlgorithm : (PairMaxDamage) how should damaged points couple
        strainAlgorithm         : (PseudoPlasticStrain) how to define the strain
        damageInCompression     : (False) optionally damage in compression as well 
                                  as tension
        criticalDamageThreshold : (4.0) prevent any nodes where Trace(D_i) exceeds 
                                  criticalDamageThreshold from setting timestep
        mask                    : (1 on all points) a field of flags: a node with 
                                  zero implies no flaws (and therefore no damage) 
                                  on that point
"""

#-------------------------------------------------------------------------------
# ProbabilisticDamageModel
#-------------------------------------------------------------------------------
ProbabilisticDamageModelGenString = """
class ProbabilisticDamageModel%(dim)s(CXXProbabilisticDamageModel%(dim)s):
    '''%(help)s'''

    def __init__(self, *args_in, **kwargs):

        args = list(args_in)

        # The material library values are in CGS, so build a units object for 
        # conversions.
        cgs = PhysicalConstants(0.01,   # unit length in m
                                0.001,  # unit mass in kg
                                1.0)    # unit time in sec

        # Arguments needed to build the damage model.
        damage_kwargs = {"nodeList"                 : None,
                         "kernel"                   : None,
                         "kWeibull"                 : None,
                         "mWeibull"                 : None,
                         "seed"                     : 48927592,
                         "minFlawsPerNode"          : 1,
                         "crackGrowthMultiplier"    : 0.4,
                         "volumeMultiplier"         : 1.0,
                         "damageCouplingAlgorithm"  : PairMaxDamage,
                         "strainAlgorithm"          : PseudoPlasticStrain,
                         "damageInCompression"      : False,
                         "criticalDamageThreshold"  : 4.0,
                         "mask"                     : None}

        # Extra arguments for our convenient constructor.
        convenient_kwargs = {"materialName"          : None,
                             "units"                 : None}

        # Check the input arguments.
        validKeys = damage_kwargs.keys() + convenient_kwargs.keys()
        for argname in kwargs:
            if not argname in validKeys:
                raise ValueError, ("ERROR: argument %%s not a valid option.\\n" %% argname +
                                   expectedUsageStringBA)

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
                raise ValueError, (("ERROR: material %%s is not in the library of material values.\\n" %% materialName) +
                                   expectedUsageStringBA)
            matprops = SpheralMaterialPropertiesLib[materialName]
            if not ("kWeibull" in matprops and "mWeibull" in matprops):
                raise ValueError, (("ERROR : material %%s does not provide the required values for kWeibull and mWeibull.\\n" %% materialName) + 
                                   expectedUsageString)
            damage_kwargs["kWeibull"] = matprops["kWeibull"]
            damage_kwargs["mWeibull"] = matprops["mWeibull"]

            # Any attempt to specify units?
            if "units" in kwargs:
                units = kwargs["units"]
                damage_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                del kwargs["units"]
            elif len(args) > 1 and isinstance(args[1], PhysicalConstants):
                units = args[1]
                damage_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                del args[1]

        # Process remaining user arguments.
        kwarg_order = ["nodeList", "kernel", "kWeibull", "mWeibull", "seed", "minFlawsPerNode", "crackGrowthMultiplier",
                       "volumeMultiplier", "damageCouplingAlgorithm", "strainAlgorithm", "damageInCompression",
                       "criticalDamageThreshold", "mask"]
        for iarg, argval in enumerate(args):
            damage_kwargs[kward_order[iarg]] = argval

        # Process any keyword arguments.  Note we already removed any deprecated keywords.
        for argname in kwargs:
            if argname in damage_kwargs:
                damage_kwargs[argname] = kwargs[argname]
            else:
                raise ValueError, (("ERROR : unknown kwarg %%s.\\n" %% argname) + expectedUsageString)

        # Build the damage model.
        CXXProbabilisticDamageModel%(dim)s.__init__(self, **damage_kwargs)
        return

"""

#-------------------------------------------------------------------------------
# Make 'em
#-------------------------------------------------------------------------------
for dim in dims:
    exec("from SpheralCompiledPackages import ProbabilisticDamageModel%id as CXXProbabilisticDamageModel%id" % (dim, dim))

    # Capture the full class help string
    save_stdout = sys.stdout
    ss = StringIO()
    sys.stdout = ss
    eval("help(CXXProbabilisticDamageModel%id)" % dim)
    sys.stdout = save_stdout
    ss.seek(0)
    class_help = ss.read()

    exec(ProbabilisticDamageModelGenString % {"dim": "%id" % dim,
                                              "help": (expectedUsageString + "\n\n Class help:\n\n" + class_help)})
