from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
SPHHydroFactoryString = """
class %(classname)s%(dim)s(SPHHydroBase%(dim)s):

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
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if WPi is None:
            WPi = W
        SPHHydroBase%(dim)s.__init__(self,
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
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "SPHHydro",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "ASPHHydro",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})
