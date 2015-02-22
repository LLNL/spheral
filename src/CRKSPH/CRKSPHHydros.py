from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(CRKSPHHydroBase%(dim)s):

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
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 momentumConserving = True):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        CRKSPHHydroBase%(dim)s.__init__(self,
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
                                      HUpdate,
                                      epsTensile,
                                      nTensile,
                                      momentumConserving)
        return
"""

#-------------------------------------------------------------------------------
# CRKSPH
#-------------------------------------------------------------------------------
exec(CRKSPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "CRKSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(CRKSPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "CRKSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(CRKSPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "CRKSPHHydro",
                              "smoothingScaleMethod" : "SPHSmoothingScale"})

#-------------------------------------------------------------------------------
# ACRKSPH
#-------------------------------------------------------------------------------
exec(CRKSPHHydroFactoryString % {"dim"                  : "1d",
                              "classname"            : "ACRKSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(CRKSPHHydroFactoryString % {"dim"                  : "2d",
                              "classname"            : "ACRKSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
exec(CRKSPHHydroFactoryString % {"dim"                  : "3d",
                              "classname"            : "ACRKSPHHydro",
                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
