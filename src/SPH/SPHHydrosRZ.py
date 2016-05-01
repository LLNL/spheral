from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()
assert 2 in dims

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
SPHHydroFactoryString = """
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
                 gradhCorrection = True,
                 XSPH = True,
                 correctVelocityGradient = True,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100)):
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
        self.zaxisBC = ReflectingBoundary2d(Plane2d(Vector2d(0.0, 0.0), Vector(0.0, 1.0)))
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
exec(SPHHydroFactoryString % {"classname"            : "SPHHydroRZ",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroFactoryString % {"classname"            : "ASPHHydroRZ",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
