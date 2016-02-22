from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.SolidSPHSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidSPHHydroFactoryString = """
class %(classname)s%(dim)s(SolidSPHHydroBase%(dim)s):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 WGrad = None,
                 filter = 0.0,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 gradhCorrection = True,
                 XSPH = True,
                 correctVelocityGradient = False,
                 sumMassDensityOverAllNodeLists = False,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if WPi is None:
            WPi = W
        if WGrad is None:
            WGrad = W
        SolidSPHHydroBase%(dim)s.__init__(self,
                                          self._smoothingScaleMethod,
                                          Q,
                                          W,
                                          WPi,
                                          WGrad,
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
    exec(SolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidSPHHydro",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidASPHHydro",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})
