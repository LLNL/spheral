import copy
from SolidSpheral import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

#-------------------------------------------------------------------------------
# Define help strings for our constructors.
#-------------------------------------------------------------------------------
expectedUsageStringBA = """
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
        strainAlgorithm     : (optional) defaults to "PsuedoPlasticStrain"
        effectiveDamageAlgorithm : (optional) defaults to "Copy".
        useDamageGradient   : (optional) defaults to "True"
        crackGrowthMultiplier : (optional) defaults to "0.4"
        flawAlgorithm       : (optional) defaults to "FullSpectrumFlaws"
        criticalDamageThreshold : (optional) defaults to 1.0
        minFlawsPerNode     : (optional) defaults to "1"
        minTotalFlaws       : (optional) defaults to "1"
"""

expectedUsageStringO = """
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
        effectiveDamageAlgorithm : (optional) defaults to "Copy".
        useDamageGradient   : (optional) defaults to "True"
        crackGrowthMultiplier : (optional) defaults to "0.4"
        flawAlgorithm       : (optional) defaults to "FullSpectrumFlaws"
        criticalDamageThreshold : (optional) defaults to 1.0
        numFlawsPerNode     : (optional) defaults to "1"
"""

#-------------------------------------------------------------------------------
# GradyKippTensorDamage : generic definition
#-------------------------------------------------------------------------------
GradyKippTensorDamageBAGenString = """
class GradyKippTensorDamageBenzAsphaug%(dim)s(TensorDamageModel%(dim)s):

    def __init__(self, *args, **kwargs):

        # Arguments needed to build the damage model.
        damage_kwargs = {"nodeList"                 : None,
                         "strainAlgorithm"          : PseudoPlasticStrain,
                         "effectiveDamageAlgorithm" : Copy,
                         "useDamageGradient"        : True,
                         "kernel"                   : None,
                         "crackGrowthMultiplier"    : 0.4,
                         "flawAlgorithm"            : FullSpectrumFlaws,
                         "criticalDamageThreshold"  : 1.0}

        # Arguments needed to build the Weibull distribution.
        weibull_kwargs = {"volume"                   : 0.0,
                          "volumeStretchFactor"      : 1.0,
                          "seed"                     : 48927595992,
                          "kWeibull"                 : None,
                          "mWeibull"                 : None,
                          "nodeList"                 : None,
                          "minFlawsPerNode"          : 1,
                          "minTotalFlaws"            : 1}

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
            assert (x in weibull_kwargs) or (x in damage_kwargs)

        # Is the user trying to use a convenient constructor?
        iarg_start = 0
        if ((len(args) > 0 and type(args[0]) == str) or
            "materialName" in kwargs):
            if len(args) > 0 and type(args[0]) == str:
                iarg_start = 1
                materialName = args[0]
            else:
                materialName = kwargs["materialName"]
            if not materialName in SpheralMaterialPropertiesLib:
                raise ValueError, (("ERROR: material %%s is not in the library of material values.\\n" %% materialName) +
                                   expectedUsageStringBA)
            matprops = SpheralMaterialPropertiesLib[materialName]
            if not ("kWeibull" in matprops and "mWeibull" in matprops):
                raise ValueError, (("ERROR : material %%s does not provide the required values for kWeibull and mWeibull.\\n" %% materialName) + 
                                   expectedUsageStringBA)
            weibull_kwargs["kWeibull"] = matprops["kWeibull"]
            weibull_kwargs["mWeibull"] = matprops["mWeibull"]

        # Process the other user arguments.
        for iarg in xrange(iarg_start, len(args)):
            argname = backCompatOrder[iarg]
            if argname in damage_kwargs:
                damage_kwargs[argname] = args[iarg]
            if argname in weibull_kwargs:
                weibull_kwargs[argname] = args[iarg]

        # Process any keyword arguments.
        for argname in kwargs:
            if argname in damage_kwargs:
                damage_kwargs[argname] = kwargs[argname]
            if argname in weibull_kwargs:
                weibull_kwargs[argname] = kwargs[argname]

        # Remember a few critical pieces of state.
        self.seed = weibull_kwargs["seed"]
        self.kWeibull = weibull_kwargs["kWeibull"]
        self.mWeibull = weibull_kwargs["mWeibull"]

        # Build the flaw distribution.
        damage_kwargs["flaws"] = weibullFlawDistributionBenzAsphaug%(dim)s(**weibull_kwargs)

        # Invoke the parent constructor.
        TensorDamageModel%(dim)s.__init__(self, **damage_kwargs)

        return

    def label(self):
        return "GradyKippTensorDamageBenzAsphaug"

    def dumpState(self,
                   file,
                   pathName):
        TensorDamageModel%(dim)s.dumpState(self, file, pathName)
        file.writeObject(self.kWeibull, pathName + "/kWeibull")
        file.writeObject(self.mWeibull, pathName + "/mWeibull")
        file.writeObject(self.seed, pathName + "/seed")
        return

    def restoreState(self,
                      file,
                      pathName):
        TensorDamageModel%(dim)s.restoreState(self, file, pathName)
        self.kWeibull = file.readObject(pathName + "/kWeibull")
        self.mWeibull = file.readObject(pathName + "/mWeibull")
        self.seed = file.readObject(pathName + "/seed")
        return

"""

#-------------------------------------------------------------------------------
# GradyKippTensorDamageOwen : generic definition
#-------------------------------------------------------------------------------
GradyKippTensorDamageOwenGenString = """
class GradyKippTensorDamageOwen%(dim)s(TensorDamageModel%(dim)s):

    def __init__(self, *args, **kwargs):

        # Arguments needed to build the damage model.
        damage_kwargs = {"nodeList"                 : None,
                         "strainAlgorithm"          : PseudoPlasticStrain,
                         "effectiveDamageAlgorithm" : Copy,
                         "useDamageGradient"        : True,
                         "kernel"                   : None,
                         "crackGrowthMultiplier"    : 0.4,
                         "flawAlgorithm"            : FullSpectrumFlaws,
                         "criticalDamageThreshold"  : 1.0}

        # Arguments needed to build the Weibull distribution.
        weibull_kwargs = {"seed"                     : 48927595992,
                          "kWeibull"                 : None,
                          "mWeibull"                 : None,
                          "nodeList"                 : None,
                          "numFlawsPerNode"          : 1,
                          "volumeMultiplier"         : 1.0}

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
                           "numFlawsPerNode"]
        for x in backCompatOrder:
            assert (x in weibull_kwargs) or (x in damage_kwargs)

        # Is the user trying to use a convenient constructor?
        iarg_start = 0
        if ((len(args) > 0 and type(args[0]) == str) or
            "materialName" in kwargs):
            if len(args) > 0 and type(args[0]) == str:
                iarg_start = 1
                materialName = args[0]
            else:
                materialName = kwargs["materialName"]
            if not materialName in SpheralMaterialPropertiesLib:
                raise ValueError, (("ERROR: material %%s is not in the library of material values.\\n" %% materialName) +
                                   expectedUsageStringO)
            matprops = SpheralMaterialPropertiesLib[materialName]
            if not ("kWeibull" in matprops and "mWeibull" in matprops):
                raise ValueError, (("ERROR : material %%s does not provide the required values for kWeibull and mWeibull.\\n" %% materialName) + 
                                   expectedUsageStringO)
            weibull_kwargs["kWeibull"] = matprops["kWeibull"]
            weibull_kwargs["mWeibull"] = matprops["mWeibull"]

        # Process the other user arguments.
        for iarg in xrange(iarg_start, len(args)):
            argname = backCompatOrder[iarg]
            if argname in damage_kwargs:
                damage_kwargs[argname] = args[iarg]
            if argname in weibull_kwargs:
                weibull_kwargs[argname] = args[iarg]

        # Process any keyword arguments.
        for argname in kwargs:
            if argname in damage_kwargs:
                damage_kwargs[argname] = kwargs[argname]
            if argname in weibull_kwargs:
                weibull_kwargs[argname] = kwargs[argname]

        # Remember a few critical pieces of state.
        self.seed = weibull_kwargs["seed"]
        self.kWeibull = weibull_kwargs["kWeibull"]
        self.mWeibull = weibull_kwargs["mWeibull"]

        # Build the flaw distribution.
        damage_kwargs["flaws"] = weibullFlawDistributionOwen%(dim)s(**weibull_kwargs)

        # Invoke the parent constructor.
        TensorDamageModel%(dim)s.__init__(self, **damage_kwargs)

        return

    def label(self):
        return "GradyKippTensorDamageOwen"

    def dumpState(self,
                   file,
                   pathName):
        TensorDamageModel%(dim)s.dumpState(self, file, pathName)
        file.writeObject(self.kWeibull, pathName + "/kWeibull")
        file.writeObject(self.mWeibull, pathName + "/mWeibull")
        file.writeObject(self.seed, pathName + "/seed")
        return

    def restoreState(self,
                      file,
                      pathName):
        TensorDamageModel%(dim)s.restoreState(self, file, pathName)
        self.kWeibull = file.readObject(pathName + "/kWeibull")
        self.mWeibull = file.readObject(pathName + "/mWeibull")
        self.seed = file.readObject(pathName + "/seed")
        return

"""

#-------------------------------------------------------------------------------
# GradyKippTensorDamageBenzAsphaug instantiations.
#-------------------------------------------------------------------------------
exec(GradyKippTensorDamageBAGenString % {"dim": "1d"})
exec(GradyKippTensorDamageBAGenString % {"dim": "2d"})
exec(GradyKippTensorDamageBAGenString % {"dim": "3d"})

#-------------------------------------------------------------------------------
# GradyKippTensorDamageOwen instantiations.
#-------------------------------------------------------------------------------
exec(GradyKippTensorDamageOwenGenString % {"dim": "1d"})
exec(GradyKippTensorDamageOwenGenString % {"dim": "2d"})
exec(GradyKippTensorDamageOwenGenString % {"dim": "3d"})

#-------------------------------------------------------------------------------
# Aliases
#-------------------------------------------------------------------------------
GradyKippTensorDamage1d = GradyKippTensorDamageBenzAsphaug1d
GradyKippTensorDamage2d = GradyKippTensorDamageBenzAsphaug2d
GradyKippTensorDamage3d = GradyKippTensorDamageBenzAsphaug3d
