from SpheralModules.Spheral.TaylorSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

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
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(TaylorSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                        "classname"            : "TaylorSPHHydro",
                                        "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(TaylorSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                        "classname"            : "TaylorASPHHydro",
                                        "smoothingScaleMethod" : "ASPHSmoothingScale"})
