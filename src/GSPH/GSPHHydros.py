from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# density-based GSPH factory string
#-------------------------------------------------------------------------------
GSPHHydroFactoryString = """
class %(classname)s%(dim)s(GSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 riemannSolver,
                 W,
                 epsDiffusionCoeff = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 gradientType = HydroAccelerationGradient,
                 densityUpdate = IntegrateDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        GSPHHydroBase%(dim)s.__init__(self,
                                     self._smoothingScaleMethod,
                                     dataBase,
                                     riemannSolver,
                                     W,
                                     epsDiffusionCoeff,
                                     cfl,
                                     useVelocityMagnitudeForDt,
                                     compatibleEnergyEvolution,
                                     evolveTotalEnergy,
                                     XSPH,
                                     correctVelocityGradient,
                                     gradientType,
                                     densityUpdate,
                                     HUpdate,
                                     epsTensile,
                                     nTensile,
                                     xmin,
                                     xmax)
        return
"""

#-------------------------------------------------------------------------------
# volume-based GSPH factory string (MFM)
#-------------------------------------------------------------------------------
MFMHydroFactoryString = """
class %(classname)s%(dim)s(MFMHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 riemannSolver,
                 W,
                 epsDiffusionCoeff = 0.0,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 gradientType = HydroAccelerationGradient,
                 densityUpdate = IntegrateDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        MFMHydroBase%(dim)s.__init__(self,
                                     self._smoothingScaleMethod,
                                     dataBase,
                                     riemannSolver,
                                     W,
                                     epsDiffusionCoeff,
                                     cfl,
                                     useVelocityMagnitudeForDt,
                                     compatibleEnergyEvolution,
                                     evolveTotalEnergy,
                                     XSPH,
                                     correctVelocityGradient,
                                     gradientType,
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
    exec(GSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "GSPHHydro",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(GSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "AGSPHHydro",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})

    exec(MFMHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "MFMHydro",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(MFMHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "AMFMHydro",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})
                                  
#-------------------------------------------------------------------------------
# GSPH convience wrapper function
#-------------------------------------------------------------------------------
def GSPH(dataBase,
        W,
        riemannSolver=None,
        specificThermalEnergyDiffusionCoefficient = 0.0,
        cfl = 0.25,
        gradientType = HydroAccelerationGradient,
        densityUpdate = IntegrateDensity,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        XSPH = False,
        correctVelocityGradient = False,
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = False,
        negativePressureInDamage = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        RZ = False):


    assert densityUpdate in (RigorousSumDensity,IntegrateDensity)

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "GSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "            If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "            which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    if ASPH:
        Constructor = eval("AGSPHHydro%id" % ndim)
    else:
        Constructor = eval("GSPHHydro%id" % ndim)

    if riemannSolver is None:
        waveSpeedMethod = eval("DavisWaveSpeed%id()" % (ndim))
        slopeLimiter = eval("VanLeerLimiter%id()" % (ndim))
        linearReconstruction=True
        riemannSolver = eval("HLLC%id(slopeLimiter,waveSpeedMethod,linearReconstruction)" % (ndim))
   
    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax

    kwargs = {"riemannSolver" : riemannSolver,
              "W" : W,
              "epsDiffusionCoeff" : specificThermalEnergyDiffusionCoefficient,
              "dataBase" : dataBase,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "gradientType" : gradientType,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}


    # Build and return the thing.
    result = Constructor(**kwargs)
    return result

#-------------------------------------------------------------------------------
# MFM convience wrapper function
#-------------------------------------------------------------------------------
def MFM(dataBase,
        W,
        riemannSolver=None,
        specificThermalEnergyDiffusionCoefficient = 0.0,
        cfl = 0.25,
        gradientType = HydroAccelerationGradient,
        densityUpdate = IntegrateDensity,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        XSPH = False,
        correctVelocityGradient = False,
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = False,
        negativePressureInDamage = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        RZ = False):

    # for now we'll just piggy back off this enum
    assert densityUpdate in (RigorousSumDensity,IntegrateDensity)

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "MFM  Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "            If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "            which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    if ASPH:
        Constructor = eval("AMFMHydro%id" % ndim)
    else:
        Constructor = eval("MFMHydro%id" % ndim)

    if riemannSolver is None:
        waveSpeedMethod = eval("DavisWaveSpeed%id()" % (ndim))
        slopeLimiter = eval("VanLeerLimiter%id()" % (ndim))
        linearReconstruction=True
        riemannSolver = eval("HLLC%id(slopeLimiter,waveSpeedMethod,linearReconstruction)" % (ndim))
   
    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax

    kwargs = {"riemannSolver" : riemannSolver,
              "W" : W,
              "epsDiffusionCoeff" : specificThermalEnergyDiffusionCoefficient,
              "dataBase" : dataBase,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "gradientType" : gradientType,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}


    # Build and return the thing.
    result = Constructor(**kwargs)
    return result

