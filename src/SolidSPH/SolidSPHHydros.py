from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.SolidSPHSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic SolidSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidSPHHydroFactoryString = """
class %(classname)s%(dim)s(SolidSPHHydroBase%(dim)s):

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
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.3,
                 nTensile = 4.0):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        SolidSPHHydroBase%(dim)s.__init__(self,
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
                                          densityUpdate,
                                          HUpdate,
                                          epsTensile,
                                          nTensile)
        return
"""

#-------------------------------------------------------------------------------
# SolidSPH
#-------------------------------------------------------------------------------
exec(SolidSPHHydroFactoryString % {"dim"                  : "1d",
                                   "classname"            : "SolidSPHHydro",
                                   "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SolidSPHHydroFactoryString % {"dim"                  : "2d",
                                   "classname"            : "SolidSPHHydro",
                                   "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SolidSPHHydroFactoryString % {"dim"                  : "3d",
                                   "classname"            : "SolidSPHHydro",
                                   "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# SolidASPH
#-------------------------------------------------------------------------------
exec(SolidSPHHydroFactoryString % {"dim"                  : "1d",
                                   "classname"            : "SolidASPHHydro",
                                   "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SolidSPHHydroFactoryString % {"dim"                  : "2d",
                                   "classname"            : "SolidASPHHydro",
                                   "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SolidSPHHydroFactoryString % {"dim"                  : "3d",
                                   "classname"            : "SolidASPHHydro",
                                   "smoothingScaleMethod" : "ASPHSmoothingScale"})
