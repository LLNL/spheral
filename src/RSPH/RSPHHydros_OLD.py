from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic RSPHHydro pattern.
#-------------------------------------------------------------------------------
RSPHHydroFactoryString = """
class %(classname)s%(dim)s(RSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 W,
                 cfl = 0.25,
                 useVelocityMagnitudeForDt = False,
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
                                      Q,
                                      W,
                                      cfl,
                                      useVelocityMagnitudeForDt,
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
    
#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SPH hydro.
#-------------------------------------------------------------------------------
def RSPH(dataBase,
         W,
         cfl = 0.25,
         useVelocityMagnitudeForDt = False,
         XSPH = True,
         correctVelocityGradient = True,
         HUpdate = IdealH,
         epsTensile = 0.0,
         nTensile = 4.0,
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
        print "SPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "           If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "           which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    # Decide on the hydro object.
    Constructor = eval("RSPHHydro%id" % ndim)

    # Dummy Artificial viscosity for now
    Cl = 1.0*(dataBase.maxKernelExtent/2.0)
    Cq = 1.0*(dataBase.maxKernelExtent/2.0)**2
    Q = eval("MonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "Q" : Q,
              "W" : W,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    #if nsolid > 0:
    #    kwargs.update({"damageRelieveRubble"      : damageRelieveRubble,
    #                   "negativePressureInDamage" : negativePressureInDamage,
    #                   "strengthInDamage"         : strengthInDamage})

    # Build and return the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    return result

