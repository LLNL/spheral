from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.BoundarySpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()
assert 2 in dims

#-------------------------------------------------------------------------------
# The generic SPHHydroRZ pattern.
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
        self.zaxisPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
        self.zaxisBC = ReflectingBoundary2d(self.zaxisPlane)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The generic SPHHydroAreaRZ pattern.
#-------------------------------------------------------------------------------
SPHHydroAreaRZFactoryString = """
class %(classname)s(SPHHydroBaseAreaRZ):

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
        SPHHydroBaseAreaRZ.__init__(self,
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
        self.zaxisPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
        self.zaxisBC = ReflectingBoundary2d(self.zaxisPlane)
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

exec(SPHHydroAreaRZFactoryString % {"classname"            : "SPHHydroAreaRZ",
                                    "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroAreaRZFactoryString % {"classname"            : "ASPHHydroAreaRZ",
                                    "smoothingScaleMethod" : "ASPHSmoothingScale"})
