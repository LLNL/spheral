from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
SPHHydroFactoryString = """
class %(classname)s%(dim)s(SPHHydroBase%(dim)s):

    def __init__(self,
                 W,
                 WPi,
                 Q,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = False,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.3,
                 nTensile = 4.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if xmin is None:
            xmin = Vector%(dim)s.one
        SPHHydroBase%(dim)s.__init__(self,
                                     self._smoothingScaleMethod,
                                     W,
                                     WPi,
                                     Q,
                                     cfl,
                                     useVelocityMagnitudeForDt,
                                     compatibleEnergyEvolution,
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
        return
"""

#-------------------------------------------------------------------------------
# SPH
#-------------------------------------------------------------------------------
exec(SPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "SPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "SPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "SPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# ASPH
#-------------------------------------------------------------------------------
exec(SPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "ASPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "ASPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "ASPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
