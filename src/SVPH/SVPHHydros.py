from SpheralModules.Spheral import *
from SpheralModules.Spheral.SVPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *

#-------------------------------------------------------------------------------
# The generic SVPHHydro pattern.
#-------------------------------------------------------------------------------
FactoryString = """
class %(classname)s%(dim)s(%(basename)sBase%(dim)s):

    def __init__(self,
                 W,
                 Q,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSVPH = True,
                 linearConsistent = False,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if xmin is None:
            xmin = Vector%(dim)s.one
        %(basename)sBase%(dim)s.__init__(self,
                                      self._smoothingScaleMethod,
                                      W,
                                      Q,
                                      cfl,
                                      useVelocityMagnitudeForDt,
                                      compatibleEnergyEvolution,
                                      XSVPH,
                                      linearConsistent,
                                      densityUpdate,
                                      HUpdate,
                                      xmin,
                                      xmax)
        return
"""

for dim in ("1d", "2d", "3d"):
    for basename in ("SVPHHydro", "SVPHFacetedHydro"):
        for prefix, smoothingScaleMethod in (("", "SPHSmoothingScale"),
                                             ("A", "ASPHSmoothingScale")):
            exec(FactoryString % {"dim"                  : dim,
                                  "basename"             : basename,
                                  "classname"            : prefix + basename,
                                  "smoothingScaleMethod" : smoothingScaleMethod})
