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
                 volumeType = CRKVoronoiVolume,
                 detectSurfaces = True,
                 detectThreshold = 0.95,
                 sweepAngle = 0.8,
                 detectRange = 0.5,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 etaMinAxis = 0.1):
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
                                   detectSurfaces,
                                   detectThreshold,
                                   sweepAngle,
                                   detectRange,
                                   epsTensile,
                                   nTensile)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
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
