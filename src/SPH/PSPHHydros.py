from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic PSPHHydro pattern.
#-------------------------------------------------------------------------------
PSPHHydroFactoryString = """
class %(classname)s%(dim)s(PSPHHydroBase%(dim)s):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
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
        if WPi is None:
            WPi = W
        PSPHHydroBase%(dim)s.__init__(self,
                                      self._smoothingScaleMethod,
                                      Q,
                                      W,
                                      WPi,
                                      filter,
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
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(PSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                   "classname"            : "PSPHHydro",
                                   "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(PSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                   "classname"            : "APSPHHydro",
                                   "smoothingScaleMethod" : "ASPHSmoothingScale"})
