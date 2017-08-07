from SpheralModules.Spheral import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.BoundarySpace import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The area-weighted SPHHydroRZ objects.
#-------------------------------------------------------------------------------
SPHHydroRZFactoryString = """
class %(classname)s(SPHHydroBaseRZ):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 sumMassDensityOverAllNodeLists = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100),
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        SPHHydroBaseRZ.__init__(self,
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
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The Garcia-Senz SPHHydroGSRZ objects.
#-------------------------------------------------------------------------------
SPHHydroGSRZFactoryString = """
class %(classname)s(SPHHydroBaseGSRZ):

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
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100),
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        SPHHydroBaseGSRZ.__init__(self,
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
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The SolidSPHHydroRZ pattern.
#-------------------------------------------------------------------------------
SolidSPHHydroRZFactoryString = """
class %(classname)s(SolidSPHHydroBaseRZ):

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
                 correctVelocityGradient = True,
                 sumMassDensityOverAllNodeLists = False,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector2d(-1e100, -1e100),
                 xmax = Vector2d( 1e100,  1e100),
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        if WGrad is None:
            WGrad = W
        SolidSPHHydroBaseRZ.__init__(self,
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
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
exec(SPHHydroRZFactoryString % {"classname"            : "SPHHydroRZ",
                                "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroRZFactoryString % {"classname"            : "ASPHHydroRZ",
                                "smoothingScaleMethod" : "ASPHSmoothingScale"})

exec(SPHHydroGSRZFactoryString % {"classname"            : "SPHHydroGSRZ",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SPHHydroGSRZFactoryString % {"classname"            : "ASPHHydroGSRZ",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})

exec(SolidSPHHydroRZFactoryString % {"classname"            : "SolidSPHHydroRZ",
                                     "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SolidSPHHydroRZFactoryString % {"classname"            : "SolidASPHHydroRZ",
                                     "smoothingScaleMethod" : "ASPHSmoothingScale"})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SPH hydro.
#-------------------------------------------------------------------------------
def SPHRZ(dataBase,
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
          xmax = ( 1e100,  1e100,  1e100),
          etaMinAxis = 0.1,
          ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    assert ndim == 2
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
            Constructor = SolidASPHHydroRZ
        else:
            Constructor = SolidSPHHydroRZ
    else:
        if ASPH:
            Constructor = ASPHHydroRZ
        else:
            Constructor = SPHHydroRZ

    # Build and return the thing.
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    return Constructor(Q = Q,
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
                       xmax = eval("Vector%id(%g, %g, %g)" % xmax),
                       etaMinAxis = etaMinAxis)

#-------------------------------------------------------------------------------
# Provide a shorthand ASPH factory.
#-------------------------------------------------------------------------------
def ASPHRZ(dataBase,
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
           xmax = ( 1e100,  1e100,  1e100),
           etaMinAxis = 0.1):
    return SPHRZ(dataBase = dataBase,
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
                 etaMinAxis = etaMinAxis,
                 ASPH = True)
