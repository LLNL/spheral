from SpheralModules.Spheral import *
from SpheralModules.Spheral.SVPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic SVPHHydro pattern.
#-------------------------------------------------------------------------------
SVPHHydroFactoryString = """
class %(classname)s%(dim)s(SVPHHydroBase%(dim)s):

    def __init__(self,
                 W,
                 Q,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSVPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if xmin is None:
            xmin = Vector%(dim)s.one
        SVPHHydroBase%(dim)s.__init__(self,
                                      self._smoothingScaleMethod,
                                      W,
                                      Q,
                                      cfl,
                                      useVelocityMagnitudeForDt,
                                      compatibleEnergyEvolution,
                                      XSVPH,
                                      densityUpdate,
                                      HUpdate,
                                      xmin,
                                      xmax)
        return
"""

#-------------------------------------------------------------------------------
# SVPH
#-------------------------------------------------------------------------------
exec(SVPHHydroFactoryString % {"dim"                  : "1d",
                               "classname"            : "SVPHHydro",
                               "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SVPHHydroFactoryString % {"dim"                  : "2d",
                               "classname"            : "SVPHHydro",
                               "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SVPHHydroFactoryString % {"dim"                  : "3d",
                               "classname"            : "SVPHHydro",
                               "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# ASVPH
#-------------------------------------------------------------------------------
exec(SVPHHydroFactoryString % {"dim"                  : "1d",
                               "classname"            : "ASVPHHydro",
                               "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SVPHHydroFactoryString % {"dim"                  : "2d",
                               "classname"            : "ASVPHHydro",
                               "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(SVPHHydroFactoryString % {"dim"                  : "3d",
                               "classname"            : "ASVPHHydro",
                               "smoothingScaleMethod" : "ASPHSmoothingScale"})
