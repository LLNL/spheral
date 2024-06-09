import sys
from SpheralCompiledPackages import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# GradyKippTensorDamageBenzAsphaug
#-------------------------------------------------------------------------------
def GKTBAFactory(ndim):
    IntField = eval("IntField{}d".format(ndim))
    TensorDamageModel = eval("TensorDamageModel{}d".format(ndim))
    weibullFlawDistributionBenzAsphaug = eval("weibullFlawDistributionBenzAsphaug{}d".format(ndim))
    class GradyKippTensorDamageBenzAsphaug(TensorDamageModel):
        def __init__(self, *args, **kwargs):
            """
GradyKippTensorDamageBenzAsphaug is constructed with the following arguments:

        materialName        : (optional) label for the material in data base to
                              lookup the Weibull parameters
        nodeList            : (required) the FluidNodeList this model should be 
                              applied to
        kWeibull            : the "k" Weibull constant -- can be looked up from 
                              materialName or provided
        mWeibull            : the "m" Weibull constant -- can be looked up from 
                              materialName or provided
        volume              : (optional) The volume to use initializing flaws.
                              Defaults to "0.0", which causes the volume to be
                              calculated from the NodeList.
        volumeStretchFactor : (optional) Multiplies per node volume, useful for
                              reduced dimensional (1D or 2D) problems.
        kernel              : (required) the interpolation kernel to use
        seed                : (optional) random number seed for flaw generation.
        strainAlgorithm     : (optional) defaults to "BenzAsphaugStrain"
        damageCouplingAlgorithm : (optional) defaults to "PairMaxDamage"
        crackGrowthMultiplier   : (optional) defaults to "0.4"
        criticalDamageThreshold : (optional) defaults to 4.0 (inactive)
        minFlawsPerNode     : (optional) defaults to "1"
        minTotalFlaws       : (optional) defaults to "1"
        mask                : (optional) a field of flags: a node with zero implies
                              do not initialize flaws on that node.  default=None
"""

            # The material library values are in CGS, so build a units object for 
            # conversions.
            cgs = PhysicalConstants(0.01,   # unit length in m
                                    0.001,  # unit mass in kg
                                    1.0)    # unit time in sec

            # Arguments needed to build the damage model.
            damage_kwargs = {"nodeList"                 : None,
                             "strainAlgorithm"          : BenzAsphaugStrain,
                             "damageCouplingAlgorithm"  : PairMaxDamage,
                             "kernel"                   : None,
                             "crackGrowthMultiplier"    : 0.4,
                             "criticalDamageThreshold"  : 4.0,
                             "damageInCompression"      : False}

            # Arguments needed to build the Weibull distribution.
            weibull_kwargs = {"volume"                   : 0.0,
                              "volumeStretchFactor"      : 1.0,
                              "seed"                     : 48927592,
                              "kWeibull"                 : None,
                              "mWeibull"                 : None,
                              "nodeList"                 : None,
                              "minFlawsPerNode"          : 1,
                              "minTotalFlaws"            : 1,
                              "mask"                     : None}

            # Extra arguments for our convenient constructor.
            convenient_kwargs = {"materialName"          : None,
                                 "units"                 : None}

            # Deprecated arguments.
            deprecated_kwargs = ["effectiveDamageAlgorithm",
                                 "useDamageGradient",
                                 "flawAlgorithm"]

            # The order of all arguments in the original constructor for these classes.
            backCompatOrder = ["nodeList",
                               "kWeibull",
                               "mWeibull",
                               "volume",
                               "volumeStretchFactor",
                               "kernel",
                               "seed",
                               "strainAlgorithm",
                               "effectiveDamageAlgorithm",
                               "useDamageGradient",
                               "crackGrowthMultiplier",
                               "flawAlgorithm",
                               "criticalDamageThreshold",
                               "minFlawsPerNode",
                               "minTotalFlaws"]
            for x in backCompatOrder:
                assert (x in weibull_kwargs) or (x in damage_kwargs) or (x in deprecated_kwargs)

            # Check the input arguments.
            validKeys = list(damage_kwargs.keys()) + list(weibull_kwargs.keys()) + list(convenient_kwargs.keys()) + deprecated_kwargs
            for argname in kwargs:
                if not argname in validKeys:
                    raise ValueError("ERROR: argument {} not a valid option.\n".format(argname))

            # Did the user try any convenient constructor operations?
            if ((len(args) > 0 and type(args[0]) == str) or
                "materialName" in kwargs):
                if len(args) > 0 and type(args[0]) == str:
                    materialName = args[0]
                    args = tuple(args[1:])                    #  del args[0]
                else:
                    materialName = kwargs["materialName"]
                    del kwargs["materialName"]
                if not materialName in SpheralMaterialPropertiesLib:
                    raise ValueError("ERROR: material {} is not in the library of material values.\n".format(materialName))
                matprops = SpheralMaterialPropertiesLib[materialName]
                if not ("kWeibull" in matprops and "mWeibull" in matprops):
                    raise ValueError("ERROR : material {} does not provide the required values for kWeibull and mWeibull.\n".format(materialName))
                weibull_kwargs["kWeibull"] = matprops["kWeibull"]
                weibull_kwargs["mWeibull"] = matprops["mWeibull"]

                # Any attempt to specify units?
                if "units" in kwargs:
                    units = kwargs["units"]
                    weibull_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                    del kwargs["units"]

            # Process the other user arguments.
            for iarg in range(len(args)):
                argname = backCompatOrder[iarg]
                if argname in deprecated_kwargs:
                    sys.stdout.write("WARNING: constructor argument {} is deprecated and will be ignored.\n".format(argname))
                if argname in damage_kwargs:
                    damage_kwargs[argname] = args[iarg]
                if argname in weibull_kwargs:
                    weibull_kwargs[argname] = args[iarg]

            # Process any keyword arguments.  Note we already removed any deprecated keywords.
            for argname in kwargs:
                if argname in deprecated_kwargs:
                    sys.stdout.write("WARNING: constructor argument {} is deprecated and will be ignored.\n".format(argname))
                if argname in damage_kwargs:
                    damage_kwargs[argname] = kwargs[argname]
                if argname in weibull_kwargs:
                    weibull_kwargs[argname] = kwargs[argname]

            # Remember a few critical pieces of state.
            self.seed = weibull_kwargs["seed"]
            self.kWeibull = weibull_kwargs["kWeibull"]
            self.mWeibull = weibull_kwargs["mWeibull"]

            # Check for any mask on generating Weibull flaws.
            if weibull_kwargs["mask"] is None:
                weibull_kwargs["mask"] = IntField("mask", damage_kwargs["nodeList"], 1)

            # Preserve the input for constructing the flaws
            self.weibull_kwargs = weibull_kwargs

            # Invoke the parent constructor.
            TensorDamageModel.__init__(self, **damage_kwargs)
            return

        ########################################################################
        def initializeProblemStartupDependencies(self,
                                                 dataBase,
                                                 state,
                                                 derivs):
            # Set the flaws
            self.weibull_kwargs["state"] = state
            self.flaws = weibullFlawDistributionBenzAsphaug(**self.weibull_kwargs)

            TensorDamageModel.initializeProblemStartupDependencies(self,
                                                                   dataBase,
                                                                   state,
                                                                   derivs)

            return

        ########################################################################
        def label(self):
            return "GradyKippTensorDamageBenzAsphaug"

    return GradyKippTensorDamageBenzAsphaug

#-------------------------------------------------------------------------------
# GradyKippTensorDamageOwen
#-------------------------------------------------------------------------------
def GKTDOFactory(ndim):
    IntField = eval("IntField{}d".format(ndim))
    TensorDamageModel = eval("TensorDamageModel{}d".format(ndim))
    weibullFlawDistributionOwen = eval("weibullFlawDistributionOwen{}d".format(ndim))
    class GradyKippTensorDamageOwen(TensorDamageModel):
        def __init__(self, *args, **kwargs):
            """
GradyKippTensorDamageOwen is constructed with the following arguments:

        materialName        : (optional) label for the material in data base to
                              lookup the Weibull parameters
        nodeList            : (required) the FluidNodeList this model should be 
                              applied to
        kWeibull            : the "k" Weibull constant -- can be looked up from 
                              materialName or provided
        mWeibull            : the "m" Weibull constant -- can be looked up from 
                              materialName or provided
        kernel              : (required) the interpolation kernel to use
        seed                : (optional) random number seed for flaw generation.
        volumeMultiplier    : (optional) Multiplies the total volume.
        strainAlgorithm     : (optional) defaults to "PsuedoPlasticStrain"
        damageCouplingAlgorithm : (optional) defaults to "PairMaxDamage"
        crackGrowthMultiplier   : (optional) defaults to "0.4"
        criticalDamageThreshold : (optional) defaults to 4.0 (inactive)
        minFlawsPerNode     : (optional) defaults to "1"
        mask                : (optional) a field of flags: a node with zero implies
                              do not initialize flaws on that node.  default=None
"""

            # The material library values are in CGS, so build a units object for 
            # conversions.
            cgs = PhysicalConstants(0.01,   # unit length in m
                                    0.001,  # unit mass in kg
                                    1.0)    # unit time in sec

            # Arguments needed to build the damage model.
            damage_kwargs = {"nodeList"                 : None,
                             "strainAlgorithm"          : PseudoPlasticStrain,
                             "damageCouplingAlgorithm"  : PairMaxDamage,
                             "kernel"                   : None,
                             "crackGrowthMultiplier"    : 0.4,
                             "criticalDamageThreshold"  : 4.0,
                             "damageInCompression"      : False}

            # Arguments needed to build the Weibull distribution.
            weibull_kwargs = {"seed"                     : 48927592,
                              "kWeibull"                 : None,
                              "mWeibull"                 : None,
                              "nodeList"                 : None,
                              "minFlawsPerNode"          : 1,
                              "volumeMultiplier"         : 1.0,
                              "mask"                     : None}

            # Extra arguments for our convenient constructor.
            convenient_kwargs = {"materialName"          : None,
                                 "units"                 : None}

            # Deprecated arguments.
            deprecated_kwargs = ["effectiveDamageAlgorithm",
                                 "useDamageGradient",
                                 "flawAlgorithm"]

            # The order of all arguments in the original constructor for these classes.
            backCompatOrder = ["nodeList",
                               "kWeibull",
                               "mWeibull",
                               "kernel",
                               "seed",
                               "volumeMultiplier",
                               "strainAlgorithm",
                               "effectiveDamageAlgorithm",
                               "useDamageGradient",
                               "crackGrowthMultiplier",
                               "flawAlgorithm",
                               "criticalDamageThreshold",
                               "minFlawsPerNode"]
            for x in backCompatOrder:
                assert (x in weibull_kwargs) or (x in damage_kwargs) or (x in deprecated_kwargs)

            # Check the input arguments.
            validKeys = list(damage_kwargs.keys()) + list(weibull_kwargs.keys()) + list(convenient_kwargs.keys()) + deprecated_kwargs
            for argname in kwargs:
                if argname in deprecated_kwargs:
                    sys.stdout.write("WARNING: constructor argument {} is deprecated and will be ignored.\n".format(argname))
                elif not argname in validKeys:
                    raise ValueError("ERROR: argument {} not a valid option.\n.".format(argname))

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
                    raise ValueError("ERROR: material {} is not in the library of material values.\n".format(materialName))
                matprops = SpheralMaterialPropertiesLib[materialName]
                if not ("kWeibull" in matprops and "mWeibull" in matprops):
                    raise ValueError("ERROR : material {} does not provide the required values for kWeibull and mWeibull.\n".format(materialName))
                weibull_kwargs["kWeibull"] = matprops["kWeibull"]
                weibull_kwargs["mWeibull"] = matprops["mWeibull"]

                # Any attempt to specify units?
                if "units" in kwargs:
                    units = kwargs["units"]
                    weibull_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                    del kwargs["units"]

            # Process the other user arguments.
            for iarg in range(len(args)):
                argname = backCompatOrder[iarg]
                if argname in deprecated_kwargs:
                    sys.stdout.write("WARNING: constructor argument {} is deprecated and will be ignored.\n".format(argname))
                if argname in damage_kwargs:
                    damage_kwargs[argname] = args[iarg]
                if argname in weibull_kwargs:
                    weibull_kwargs[argname] = args[iarg]

            # Process any keyword arguments.
            for argname in kwargs:
                if argname in deprecated_kwargs:
                    sys.stdout.write("WARNING: constructor argument {} is deprecated and will be ignored.\n".format(argname))
                if argname in damage_kwargs:
                    damage_kwargs[argname] = kwargs[argname]
                if argname in weibull_kwargs:
                    weibull_kwargs[argname] = kwargs[argname]

            # Remember a few critical pieces of state.
            self.seed = weibull_kwargs["seed"]
            self.kWeibull = weibull_kwargs["kWeibull"]
            self.mWeibull = weibull_kwargs["mWeibull"]

            # Check for any mask on generating Weibull flaws.
            if weibull_kwargs["mask"] is None:
                weibull_kwargs["mask"] = IntField("mask", damage_kwargs["nodeList"], 1)

            # Preserve the input for constructing the flaws, and build dummy flaws for now
            self.weibull_kwargs = weibull_kwargs

            # Invoke the parent constructor.
            TensorDamageModel.__init__(self, **damage_kwargs)

            return

        ########################################################################
        def initializeProblemStartupDependencies(self,
                                                 dataBase,
                                                 state,
                                                 derivs):
            # Set the flaws
            self.weibull_kwargs["state"] = state
            self.flaws = weibullFlawDistributionOwen(**self.weibull_kwargs)

            TensorDamageModel.initializeProblemStartupDependencies(self,
                                                                   dataBase,
                                                                   state,
                                                                   derivs)

            return

        ########################################################################
        def label(self):
            return "GradyKippTensorDamageOwen"

    return GradyKippTensorDamageOwen

#-------------------------------------------------------------------------------
# Make 'em
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("GradyKippTensorDamageBenzAsphaug{ndim}d = GKTBAFactory({ndim})".format(ndim=ndim))
    exec("GradyKippTensorDamageOwen{ndim}d = GKTDOFactory({ndim})".format(ndim=ndim))

    # Aliases
    exec("GradyKippTensorDamage{ndim}d = GradyKippTensorDamageBenzAsphaug{ndim}d".format(ndim=ndim))
