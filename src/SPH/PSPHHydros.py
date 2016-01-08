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
                 XSPH = True,
                 correctVelocityGradient = True,
                 evolveTotalEnergy = False,
                 HopkinsConductivity = False,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
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
                                      XSPH,
                                      correctVelocityGradient,
                                      evolveTotalEnergy,
                                      HopkinsConductivity,
                                      sumMassDensityOverAllNodeLists,
                                      densityUpdate,
                                      HUpdate,
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
