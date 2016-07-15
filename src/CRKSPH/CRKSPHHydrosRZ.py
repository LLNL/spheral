from SpheralModules.Spheral import *
from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.KernelSpace import *
from SpheralModules.Spheral.BoundarySpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroRZFactoryString = """
class %(classname)s(CRKSPHHydroBaseRZ):

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
                 volumeType = CRKSumVolume,
                 epsTensile = 0.0,
                 nTensile = 4.0):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        CRKSPHHydroBaseRZ.__init__(self,
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
                                   epsTensile,
                                   nTensile)
        self.zaxisPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
        self.zaxisBC = ReflectingBoundary2d(self.zaxisPlane)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(CRKSPHHydroRZFactoryString % {"classname"            : "CRKSPHHydroRZ",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(CRKSPHHydroRZFactoryString % {"classname"            : "ACRKSPHHydroRZ",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})
