from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
RSPHHydroFactoryString = """
class %(classname)s%(dim)s(RSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 W,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        RSPHHydroBase%(dim)s.__init__(self,
                                     self._smoothingScaleMethod,
                                     dataBase,
                                     W,
                                     cfl,
                                     useVelocityMagnitudeForDt,
                                     compatibleEnergyEvolution,
                                     evolveTotalEnergy,
                                     XSPH,
                                     correctVelocityGradient,
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
    exec(RSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "RSPHHydro",
                                  "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(RSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                  "classname"            : "ARSPHHydro",
                                  "smoothingScaleMethod" : "ASPHSmoothingScale"})


#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SPH hydro.
#-------------------------------------------------------------------------------
def RSPH(dataBase,
        W,
        cfl = 0.25,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        XSPH = True,
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

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "RSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "            If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "            which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    if ASPH:
        Constructor = eval("ARSPHHydro%id" % ndim)
    else:
        Constructor = eval("RSPHHydro%id" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"W" : W,
              "dataBase" : dataBase,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}


    # Build and return the thing.
    result = Constructor(**kwargs)
    return result
