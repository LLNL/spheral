from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidSPHHydro pattern.
#-------------------------------------------------------------------------------
FSISolidSPHHydroFactoryString = """
class %(classname)s%(dim)s(FSISolidSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 W,
                 WPi = None,
                 WGrad = None,
                 alpha = 1.00,
                 diffusionCoefficient = 0.0,        
                 interfaceMethod = 0,         
                 decoupledNodeLists = None,
                 sumDensityNodeLists = None,
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
                 damageRelieveRubble = False,
                 negativePressureInDamage = False,
                 strengthInDamage = False,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        
        WPi = W
        WGrad = W

        FSISolidSPHHydroBase%(dim)s.__init__(self,
                                          self._smoothingScaleMethod,
                                          dataBase,
                                          Q,
                                          W,
                                          WPi,
                                          WGrad,
                                          alpha,
                                          diffusionCoefficient,        
                                          interfaceMethod,         
                                          decoupledNodeLists,
                                          sumDensityNodeLists,
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
                                          damageRelieveRubble,
                                          negativePressureInDamage,
                                          strengthInDamage,
                                          xmin,
                                          xmax)
        return
"""

for dim in dims:
    #exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
    #                              "classname"            : "FSISPHHydro",
    #                              "smoothingScaleMethod" : "SPHSmoothingScale"})
    #exec(SPHHydroFactoryString % {"dim"                  : "%id" % dim,
    #                              "classname"            : "AFSISPHHydro",
    #                              "smoothingScaleMethod" : "ASPHSmoothingScale"})
    
    exec(FSISolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "FSISolidSPHHydro",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(FSISolidSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "FSISolidASPHHydro",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})


def FSISPH(dataBase,
        W,
        WPi = None,
        Q = None,
        filter = 0.0,
        cfl = 0.25,
        alpha = 1.00,                 
        diffusionCoefficient=0.0,        
        interfaceMethod=0,         
        decoupledNodeLists=None,
        sumDensityNodeLists=None,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        gradhCorrection = False,
        XSPH = False,
        correctVelocityGradient = True,
        sumMassDensityOverAllNodeLists = True,
        densityUpdate = IntegrateDensity,
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

    # default to fully coupled interfaces 
    if decoupledNodeLists is None:
        decoupledNodeLists = vector_of_int([0]*dataBase.numNodeLists)
    elif isinstance(decoupledNodeLists,list):
        if len(decoupledNodeLists) == dataBase.numFluidNodeLists:
            assert (isinstance(decoupledNodeLists[0],bool) or isinstance(decoupledNodeLists[0],int))
            decoupledNodeLists = vector_of_int(decoupledNodeLists)
        else:
            raise RuntimeError, "decoupledNodeLists - must be list of length equal to numNodeLists."
    elif not isinstance(sumDensityNodeList,vector_of_int):
        raise RuntimeError, "decoupledNodeLists - must be list of int 1/0 or bool."

    # if we're using the bulk modulus interface method its going to override 
    # and user defined decoupledNodeLists to fully coupled
    if interfaceMethod == 0 and any(decoupledNodeLists) == 1:
        print "interfaceMethod = 0 requires full coupling of interface nodes. decoupledNodeLists are going to be ignored."
        decoupledNodeLists *= 0

    # default to integrate density
    if sumDensityNodeLists is None:
        sumDensityNodeLists = vector_of_int([0]*dataBase.numNodeLists)
    elif isinstance(sumDensityNodeLists,list):
        if len(sumDensityNodeLists) == dataBase.numFluidNodeLists:
            assert (isinstance(sumDensityNodeLists[0],bool) or isinstance(sumDensityNodeLists[0],int))
            sumDensityNodeLists = vector_of_int(sumDensityNodeLists)
        else:
            raise RuntimeError, "sumDensityNodeLists - must be list of length equal to numNodeLists."
    elif not isinstance(sumDensityNodeList,vector_of_int):
        raise RuntimeError, "sumDensityNodeLists - must be list of int 1/0 or bool."

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
    if RZ:

        print "SPH Error: RZ isn't set up for the FSI implementations yet"
        raise RuntimeError, "FZ inavailable for FSI."

        # RZ ----------------------------------------
        #if nsolid > 0:
        #    if ASPH:
        #        Constructor = SolidASPHHydroRZ
        #    else:
        #        Constructor = SolidSPHHydroRZ
        #else:
        #    if ASPH:
        #        Constructor = ASPHHydroRZ
        #    else:
        #        Constructor = SPHHydroRZ

    else:
        # Cartesian ---------------------------------
        if nsolid > 0:
            if ASPH:
                Constructor = eval("FSISolidASPHHydro%id" % ndim)
            else:
                Constructor = eval("FSISolidSPHHydro%id" % ndim)
        else:
            raise RuntimeError, "fluid version not available."
        #    if ASPH:
        #        Constructor = eval("ASPHHydro%id" % ndim)
        #    else:
        #        Constructor = eval("SPHHydro%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 1.0*(dataBase.maxKernelExtent/2.0)
        Cq = 1.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("MonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"W" : W,
              "WPi" : WPi,
              "dataBase" : dataBase,
              "Q" : Q,
              "filter" : filter,
              "cfl" : cfl,
              "alpha" : alpha,
              "diffusionCoefficient" : diffusionCoefficient,        
              "interfaceMethod" : interfaceMethod,         
              "decoupledNodeLists" : decoupledNodeLists,
              "sumDensityNodeLists" : sumDensityNodeLists,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "gradhCorrection" : gradhCorrection,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "sumMassDensityOverAllNodeLists" : sumMassDensityOverAllNodeLists,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    if nsolid > 0:
        kwargs.update({"damageRelieveRubble"      : damageRelieveRubble,
                       "negativePressureInDamage" : negativePressureInDamage,
                       "strengthInDamage"         : strengthInDamage})

    # Build and return the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    return result

#-------------------------------------------------------------------------------
# Provide a shorthand ASPH factory.
#-------------------------------------------------------------------------------
def AFSISPH(dataBase,
         W,
         WPi = None,
         Q = None,
         filter = 0.0,
         cfl = 0.25,
         alpha=1.00,
         diffusionCoefficient=0.0,        
         interfaceMethod=0,         
         decoupledNodeLists=None,
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
         damageRelieveRubble = False,
         negativePressureInDamage = False,
         strengthInDamage = False,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100)):
    return FSISPH(dataBase = dataBase,
               W = W,
               WPi = WPi,
               Q = Q,
               filter = filter,
               cfl = cfl,
               alpha = alpha,
               diffusionCoefficient= diffusionCoefficient,        
               interfaceMethod = interfaceMethod,         
               decoupledNodeLists = decoupledNodeLists,
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
               damageRelieveRubble = damageRelieveRubble,
               negativePressureInDamage = negativePressureInDamage,
               strengthInDamage = strengthInDamage,
               xmin = xmin,
               xmax = xmax,
               ASPH = True)
