from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.KernelSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(CRKSPHHydroBase%(dim)s):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 correctionOrder = LinearOrder,
                 volumeType = SumVolume,
                 epsTensile = 0.0,
                 nTensile = 4.0):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if WPi is None:
            WPi = W
        CRKSPHHydroBase%(dim)s.__init__(self,
                                        self._smoothingScaleMethod,
                                        Q,
                                        W,
                                        WPi,
                                        filter,
                                        cfl,
                                        useVelocityMagnitudeForDt,
                                        compatibleEnergyEvolution,
                                        XSPH,
                                        densityUpdate,
                                        HUpdate,
                                        correctionOrder,
                                        volumeType,
                                        epsTensile,
                                        nTensile)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "CRKSPHHydro",
                                     "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "ACRKSPHHydro",
                                     "smoothingScaleMethod" : "ASPHSmoothingScale"})
