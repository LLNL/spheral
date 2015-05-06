from SpheralModules.Spheral.TaylorSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic TaylorSPHHydro pattern.
#-------------------------------------------------------------------------------
TaylorSPHHydroFactoryString = """
class %(classname)s%(dim)s(TaylorSPHHydroBase%(dim)s):

    def __init__(self,
                 W,
                 Q,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSPH = True,
                 HUpdate = IdealH):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        TaylorSPHHydroBase%(dim)s.__init__(self,
                                           self._smoothingScaleMethod,
                                           W,
                                           Q,
                                           cfl,
                                           useVelocityMagnitudeForDt,
                                           compatibleEnergyEvolution,
                                           XSPH,
                                           HUpdate)
        return
"""

#-------------------------------------------------------------------------------
# TaylorSPH
#-------------------------------------------------------------------------------
exec(TaylorSPHHydroFactoryString % {"dim"                  : "1d",
                                    "classname"            : "TaylorSPHHydro",
                                    "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(TaylorSPHHydroFactoryString % {"dim"                  : "2d",
                                    "classname"            : "TaylorSPHHydro",
                                    "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(TaylorSPHHydroFactoryString % {"dim"                  : "3d",
                                    "classname"            : "TaylorSPHHydro",
                                    "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# TaylorASPH
#-------------------------------------------------------------------------------
exec(TaylorSPHHydroFactoryString % {"dim"                  : "1d",
                                    "classname"            : "TaylorASPHHydro",
                                    "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(TaylorSPHHydroFactoryString % {"dim"                  : "2d",
                                    "classname"            : "TaylorASPHHydro",
                                    "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(TaylorSPHHydroFactoryString % {"dim"                  : "3d",
                                    "classname"            : "TaylorASPHHydro",
                                    "smoothingScaleMethod" : "ASPHSmoothingScale"})
