from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.KernelSpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(CRKSPHHydroBase%(dim)s):

    def __init__(self,
                 W,
                 WPi,
                 Q,
                 Wfilter = TableKernel%(dim)s(NBSplineKernel%(dim)s(7),1000),
                 filter = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        CRKSPHHydroBase%(dim)s.__init__(self,
                                        self._smoothingScaleMethod,
                                        W,
                                        WPi,
                                        Q,
                                        Wfilter,
                                        filter,
                                        cfl,
                                        useVelocityMagnitudeForDt,
                                        compatibleEnergyEvolution,
                                        XSPH,
                                        densityUpdate,
                                        HUpdate,
                                        epsTensile,
                                        nTensile)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "CRKSPHHydro",
                                     "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "ACRKSPHHydro",
                                     "smoothingScaleMethod" : "ASPHSmoothingScale"})
