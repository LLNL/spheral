from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.KernelSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidCRKSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidCRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(SolidCRKSPHHydroBase%(dim)s):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 correctionOrder = LinearOrder,
                 volumeType = CRKVoronoiVolume,
                 detectSurfaces = True,
                 detectThreshold = 0.05,
                 sweepAngle = 0.8,
                 detectRange = 1.0,
                 epsTensile = 0.0,
                 nTensile = 4.0):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if WPi is None:
            WPi = W
        SolidCRKSPHHydroBase%(dim)s.__init__(self,
                                             self._smoothingScaleMethod,
                                             Q,
                                             W,
                                             WPi,
                                             filter,
                                             cfl,
                                             useVelocityMagnitudeForDt,
                                             compatibleEnergyEvolution,
                                             evolveTotalEnergy,
                                             XSPH,
                                             densityUpdate,
                                             HUpdate,
                                             correctionOrder,
                                             volumeType,
                                             detectSurfaces,
                                             detectThreshold,
                                             sweepAngle,
                                             detectRange,
                                             epsTensile,
                                             nTensile)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(SolidCRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                          "classname"            : "SolidCRKSPHHydro",
                                          "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidCRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                          "classname"            : "SolidACRKSPHHydro",
                                          "smoothingScaleMethod" : "ASPHSmoothingScale"})
