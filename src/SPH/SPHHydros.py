from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.ArtificialViscositySpace import *
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
                 cfl = 0.25,
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
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = True,
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
    exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "SPHHydro",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "ASPHHydro",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})

    exec(SolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidSPHHydro",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidASPHHydro",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SPH hydro.
#-------------------------------------------------------------------------------
def SPH(dataBase,
        W,
        WPi = None,
        Q = None,
        filter = 0.0,
        cfl = 0.25,
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
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "SPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "           If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "           which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    # Decide on the hydro object.
    if nsolid > 0:
        if ASPH:
            Constructor = eval("SolidASPHHydro%id" % ndim)
        else:
            Constructor = eval("SolidSPHHydro%id" % ndim)
    else:
        if ASPH:
            Constructor = eval("ASPHHydro%id" % ndim)
        else:
            Constructor = eval("SPHHydro%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Q = eval("MonaghanGingoldViscosity%id()" % ndim)

    # Build and return the thing.
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    result = Constructor(Q = Q,
                         W = W,
                         WPi = WPi,
                         filter = filter,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         evolveTotalEnergy = evolveTotalEnergy,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         correctVelocityGradient = correctVelocityGradient,
                         sumMassDensityOverAllNodeLists = sumMassDensityOverAllNodeLists,
                         densityUpdate = densityUpdate,
                         HUpdate = HUpdate,
                         epsTensile = epsTensile,
                         nTensile = nTensile,
                         xmin = eval("Vector%id(%g, %g, %g)" % xmin),
                         xmax = eval("Vector%id(%g, %g, %g)" % xmax))
    result.Q = Q
    return result

#-------------------------------------------------------------------------------
# Provide a shorthand ASPH factory.
#-------------------------------------------------------------------------------
def ASPH(dataBase,
         Q,
         W,
         WPi = None,
         filter = 0.0,
         cfl = 0.25,
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
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100)):
    return SPH(dataBase = dataBase,
               Q = Q,
               W = W,
               WPi = WPi,
               filter = filter,
               cfl = cfl,
               useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
               compatibleEnergyEvolution = compatibleEnergyEvolution,
               evolveTotalEnergy = evolveTotalEnergy,
               gradhCorrection = gradhCorrection,
               XSPH = XSPH,
               correctVelocityGradient = correctVelocityGradient,
               sumMassDensityOverAllNodeLists = sumMassDensityOverAllNodeLists,
               densityUpdate = densityUpdate,
               HUpdate = HUpdate,
               epsTensile = epsTensile,
               nTensile = nTensile,
               xmin = xmin,
               xmax = xmax,
               ASPH = True)
