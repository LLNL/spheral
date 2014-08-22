from SpheralModules.Spheral.CSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic CSPHHydro pattern.
#-------------------------------------------------------------------------------
CSPHHydroFactoryString = """
class %(classname)s%(dim)s(CSPHHydroBase%(dim)s):

    def __init__(self,
                 W,
                 WPi,
                 Q,
                 filter = 0.1,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        CSPHHydroBase%(dim)s.__init__(self,
                                      self._smoothingScaleMethod,
                                      W,
                                      WPi,
                                      Q,
                                      filter,
                                      cfl,
                                      useVelocityMagnitudeForDt,
                                      compatibleEnergyEvolution,
                                      XSPH,
                                      densityUpdate,
                                      HUpdate)
        return
"""

#-------------------------------------------------------------------------------
# CSPH
#-------------------------------------------------------------------------------
exec(CSPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "CSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(CSPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "CSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(CSPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "CSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# ACSPH
#-------------------------------------------------------------------------------
exec(CSPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "ACSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(CSPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "ACSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(CSPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "ACSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
