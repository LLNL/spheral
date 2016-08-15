from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.BoundarySpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The area-weighted SPHHydroRZ objects.
#-------------------------------------------------------------------------------
SPHHydroRZFactoryString = """
class %(classname)s(SPHHydroBaseRZ):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100),
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        SPHHydroBaseRZ.__init__(self,
                                self._smoothingScaleMethod,
                                Q,
                                W,
                                WPi,
                                filter,
                                cfl,
                                useVelocityMagnitudeForDt,
                                compatibleEnergyEvolution,
                                evolveTotalEnergy,
                                gradhCorrection,
                                XSPH,
                                correctVelocityGradient,
                                sumMassDensityOverAllNodeLists,
                                densityUpdate,
                                HUpdate,
                                epsTensile,
                                nTensile,
                                xmin,
                                xmax)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The Garcia-Senz SPHHydroGSRZ objects.
#-------------------------------------------------------------------------------
SPHHydroGSRZFactoryString = """
class %(classname)s(SPHHydroBaseGSRZ):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 gradhCorrection = True,
                 XSPH = True,
                 correctVelocityGradient = True,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100),
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        SPHHydroBaseGSRZ.__init__(self,
                                  self._smoothingScaleMethod,
                                  Q,
                                  W,
                                  WPi,
                                  filter,
                                  cfl,
                                  useVelocityMagnitudeForDt,
                                  compatibleEnergyEvolution,
                                  evolveTotalEnergy,
                                  gradhCorrection,
                                  XSPH,
                                  correctVelocityGradient,
                                  sumMassDensityOverAllNodeLists,
                                  densityUpdate,
                                  HUpdate,
                                  epsTensile,
                                  nTensile,
                                  xmin,
                                  xmax)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
exec(SPHHydroRZFactoryString % {"classname"            : "SPHHydroRZ",
                                "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroRZFactoryString % {"classname"            : "ASPHHydroRZ",
                                "smoothingScaleMethod" : "ASPHSmoothingScale"})

exec(SPHHydroGSRZFactoryString % {"classname"            : "SPHHydroGSRZ",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroGSRZFactoryString % {"classname"            : "ASPHHydroGSRZ",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})
