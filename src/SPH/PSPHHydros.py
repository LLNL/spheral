from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic PSPHHydro pattern.
#-------------------------------------------------------------------------------
PSPHHydroFactoryString = """
class %(classname)s%(dim)s(PSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 W,
                 WPi = None,
                 filter = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 correctVelocityGradient = True,
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
                                      dataBase,
                                      Q,
                                      W,
                                      WPi,
                                      filter,
                                      cfl,
                                      useVelocityMagnitudeForDt,
                                      compatibleEnergyEvolution,
                                      evolveTotalEnergy,
                                      XSPH,
                                      correctVelocityGradient,
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

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate PSPH hydro.
#-------------------------------------------------------------------------------
def PSPH(dataBase,
         W,
         WPi = None,
         Q = None,
         filter = 0.0,
         cfl = 0.25,
         useVelocityMagnitudeForDt = False,
         compatibleEnergyEvolution = True,
         evolveTotalEnergy = False,
         XSPH = True,
         correctVelocityGradient = True,
         sumMassDensityOverAllNodeLists = True,
         densityUpdate = RigorousSumDensity,
         HUpdate = IdealH,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100),
         ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "PSPH Warning: you have provided solid NodeLists, but PSPH currently does not have a solid option."
        print "              The fluid limit will be provided for now."

    # Decide on the hydro object.
    if ASPH:
        Constructor = eval("APSPHHydro%id" % ndim)
    else:
        Constructor = eval("PSPHHydro%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 2.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Build and return the thing.
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    result =  Constructor(dataBase = dataBase,
                          Q = Q,
                          W = W,
                          WPi = WPi,
                          filter = filter,
                          cfl = cfl,
                          useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                          compatibleEnergyEvolution = compatibleEnergyEvolution,
                          evolveTotalEnergy = evolveTotalEnergy,
                          XSPH = XSPH,
                          correctVelocityGradient = correctVelocityGradient,
                          sumMassDensityOverAllNodeLists = sumMassDensityOverAllNodeLists,
                          densityUpdate = densityUpdate,
                          HUpdate = HUpdate,
                          xmin = eval("Vector%id(%g, %g, %g)" % xmin),
                          xmax = eval("Vector%id(%g, %g, %g)" % xmax))
    result.Q = Q
    return result

#-------------------------------------------------------------------------------
# Provide a shorthand APSPH factory.
#-------------------------------------------------------------------------------
def PASPH(dataBase,
          W,
          WPi = None,
          Q = None,
          filter = 0.0,
          cfl = 0.25,
          useVelocityMagnitudeForDt = False,
          compatibleEnergyEvolution = True,
          evolveTotalEnergy = False,
          XSPH = True,
          correctVelocityGradient = True,
          sumMassDensityOverAllNodeLists = True,
          densityUpdate = RigorousSumDensity,
          HUpdate = IdealH,
          xmin = (-1e100, -1e100, -1e100),
          xmax = ( 1e100,  1e100,  1e100)):
    return PSPH(dataBase = dataBase,
                W = W,
                WPi = WPi,
                Q = Q,
                filter = filter,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                sumMassDensityOverAllNodeLists = sumMassDensityOverAllNodeLists,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                xmin = xmin,
                xmax = xmax,
                ASPH = True)
