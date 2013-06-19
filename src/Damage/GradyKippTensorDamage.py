from SolidSpheral import *

#-------------------------------------------------------------------------------
# GradyKippTensorDamage : generic definition
#-------------------------------------------------------------------------------
GradyKippTensorDamageBAGenString = """
class GradyKippTensorDamageBenzAsphaug%(dim)s(TensorDamageModel%(dim)s):

    def __init__(self,
                 nodeList,
                 kWeibull,
                 mWeibull,
                 volume,
                 volumeStretchFactor,
                 kernel,
                 seed,
                 strainAlgorithm = PseudoPlasticStrain,
                 effectiveDamageAlgorithm = Copy,
                 useDamageGradient = True,
                 crackGrowthMultiplier = 0.4,
                 flawAlgorithm = SampledFlaws,
                 criticalDamageThreshold = 1.0,
                 minFlawsPerNode = 1,
                 minTotalFlaws = 1):
        self.seed = seed
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        TensorDamageModel%(dim)s.__init__(self,
                                          nodeList,
                                          strainAlgorithm,
                                          effectiveDamageAlgorithm,
                                          useDamageGradient,
                                          kernel,
                                          crackGrowthMultiplier,
                                          flawAlgorithm,
                                          criticalDamageThreshold,
                                          weibullFlawDistributionBenzAsphaug%(dim)s(volume,
                                                                                    volumeStretchFactor,
                                                                                    seed,
                                                                                    kWeibull,
                                                                                    mWeibull,
                                                                                    nodeList,
                                                                                    minFlawsPerNode,
                                                                                    minTotalFlaws))
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

    def __init__(self,
                 nodeList,
                 kWeibull,
                 mWeibull,
                 kernel,
                 seed,
                 volumeMultiplier = 1.0,
                 strainAlgorithm = PseudoPlasticStrain,
                 effectiveDamageAlgorithm = Copy,
                 useDamageGradient = True,
                 crackGrowthMultiplier = 0.4,
                 flawAlgorithm = SampledFlaws,
                 criticalDamageThreshold = 1.0,
                 numFlawsPerNode = 1):
        self.seed = seed
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        TensorDamageModel%(dim)s.__init__(self,
                                          nodeList,
                                          strainAlgorithm,
                                          effectiveDamageAlgorithm,
                                          useDamageGradient,
                                          kernel,
                                          crackGrowthMultiplier,
                                          flawAlgorithm,
                                          criticalDamageThreshold,
                                          weibullFlawDistributionOwen%(dim)s(seed,
                                                                             kWeibull,
                                                                             mWeibull,
                                                                             nodeList,
                                                                             numFlawsPerNode,
                                                                             volumeMultiplier))
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
