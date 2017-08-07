from SpheralModules.Spheral import *
from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.KernelSpace import *
from SpheralModules.Spheral.BoundarySpace import *

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroRZFactoryString = """
class %(classname)s(CRKSPHHydroBaseRZ):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 1.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 correctionOrder = LinearOrder,
                 volumeType = CRKVoronoiVolume,
                 detectSurfaces = False,
                 detectThreshold = 0.05,
                 sweepAngle = 0.8,
                 detectRange = 1.0,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        CRKSPHHydroBaseRZ.__init__(self,
                                   self._smoothingScaleMethod,
                                   Q,
                                   W,
                                   WPi,
                                   filter,
                                   cfl,
                                   useVelocityMagnitudeForDt,
                                   compatibleEnergyEvolution,
                                   evolveTotalEnergy,
                                   XSPH,
                                   densityUpdate,
                                   HUpdate,
                                   correctionOrder,
                                   volumeType,
                                   detectSurfaces,
                                   detectThreshold,
                                   sweepAngle,
                                   detectRange,
                                   epsTensile,
                                   nTensile)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The generic SolidCRKSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidCRKSPHHydroRZFactoryString = """
class %(classname)s(SolidCRKSPHHydroBaseRZ):

    def __init__(self,
                 Q,
                 W,
                 WPi = None,
                 filter = 1.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 correctionOrder = LinearOrder,
                 volumeType = CRKVoronoiVolume,
                 detectSurfaces = False,
                 detectThreshold = 0.05,
                 sweepAngle = 0.8,
                 detectRange = 1.0,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 etaMinAxis = 0.1):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        if WPi is None:
            WPi = W
        SolidCRKSPHHydroBaseRZ.__init__(self,
                                        self._smoothingScaleMethod,
                                        Q,
                                        W,
                                        WPi,
                                        filter,
                                        cfl,
                                        useVelocityMagnitudeForDt,
                                        compatibleEnergyEvolution,
                                        evolveTotalEnergy,
                                        XSPH,
                                        densityUpdate,
                                        HUpdate,
                                        correctionOrder,
                                        volumeType,
                                        detectSurfaces,
                                        detectThreshold,
                                        sweepAngle,
                                        detectRange,
                                        epsTensile,
                                        nTensile)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
exec(CRKSPHHydroRZFactoryString % {"classname"            : "CRKSPHHydroRZ",
                                   "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(CRKSPHHydroRZFactoryString % {"classname"            : "ACRKSPHHydroRZ",
                                   "smoothingScaleMethod" : "ASPHSmoothingScale"})

exec(SolidCRKSPHHydroRZFactoryString % {"classname"            : "SolidCRKSPHHydroRZ",
                                        "smoothingScaleMethod" : "SPHSmoothingScale"})
exec(SolidCRKSPHHydroRZFactoryString % {"classname"            : "SolidACRKSPHHydroRZ",
                                        "smoothingScaleMethod" : "ASPHSmoothingScale"})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate CRKSPH hydro.
#-------------------------------------------------------------------------------
def CRKSPHRZ(dataBase,
             Q,
             W,
             WPi = None,
             filter = 0.0,
             cfl = 0.25,
             useVelocityMagnitudeForDt = False,
             compatibleEnergyEvolution = True,
             evolveTotalEnergy = False,
             XSPH = True,
             densityUpdate = RigorousSumDensity,
             HUpdate = IdealH,
             correctionOrder = LinearOrder,
             volumeType = CRKVoronoiVolume,
             detectSurfaces = False,
             detectThreshold = 0.05,
             sweepAngle = 0.8,
             detectRange = 1.0,
             epsTensile = 0.0,
             nTensile = 4.0,
             ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    assert ndim == 2
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "CRKSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "             If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "             which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    # Decide on the hydro object.
    if nsolid > 0:
        if ASPH:
            Constructor = SolidACRKSPHHydroRZ
        else:
            Constructor = SolidCRKSPHHydroRZ
    else:
        if ASPH:
            Constructor = ACRKSPHHydroRZ
        else:
            Constructor = CRKSPHHydroRZ

    # Build and return the thing.
    return Constructor(Q = Q,
                       W = W,
                       WPi = WPi,
                       filter = filter,
                       cfl = cfl,
                       useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                       compatibleEnergyEvolution = compatibleEnergyEvolution,
                       evolveTotalEnergy = evolveTotalEnergy,
                       XSPH = XSPH,
                       densityUpdate = densityUpdate,
                       HUpdate = HUpdate,
                       correctionOrder = correctionOrder,
                       volumeType = volumeType,
                       detectSurfaces = detectSurfaces,
                       detectThreshold = detectThreshold,
                       sweepAngle = sweepAngle,
                       detectRange = detectRange,
                       epsTensile = epsTensile,
                       nTensile = nTensile)

#-------------------------------------------------------------------------------
# Provide a shorthand ACRKSPH factory.
#-------------------------------------------------------------------------------
def ACRKSPHRZ(dataBase,
              Q,
              W,
              WPi = None,
              filter = 0.0,
              cfl = 0.25,
              useVelocityMagnitudeForDt = False,
              compatibleEnergyEvolution = True,
              evolveTotalEnergy = False,
              XSPH = True,
              densityUpdate = RigorousSumDensity,
              HUpdate = IdealH,
              correctionOrder = LinearOrder,
              volumeType = CRKVoronoiVolume,
              detectSurfaces = False,
              detectThreshold = 0.05,
              sweepAngle = 0.8,
              detectRange = 1.0,
              epsTensile = 0.0,
              nTensile = 4.0):
    return CRKSPHRZ(dataBase = dataBase,
                    Q = Q,
                    W = W,
                    WPi = WPi,
                    filter = filter,
                    cfl = cfl,
                    useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                    compatibleEnergyEvolution = compatibleEnergyEvolution,
                    evolveTotalEnergy = evolveTotalEnergy,
                    XSPH = XSPH,
                    densityUpdate = densityUpdate,
                    HUpdate = HUpdate,
                    correctionOrder = correctionOrder,
                    volumeType = volumeType,
                    detectSurfaces = detectSurfaces,
                    detectThreshold = detectThreshold,
                    sweepAngle = sweepAngle,
                    detectRange = detectRange,
                    epsTensile = epsTensile,
                    nTensile = nTensile,
                    ASPH = True)
