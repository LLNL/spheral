from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic solidFSISPHHydro pattern.
#-------------------------------------------------------------------------------
SolidFSISPHHydroFactoryString = """
class %(classname)s%(dim)s(SolidFSISPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 slides,
                 W,
                 filter = 0.0,
                 cfl = 0.5,
                 surfaceForceCoefficient = 0.0,
                 densityStabilizationCoefficient = 0.0, 
                 specificThermalEnergyDiffusionCoefficient = 0.0, 
                 xsphCoefficient=0.0,
                 interfaceMethod=HLLCInterface,  
                 kernelAveragingMethod=NeverAverageKernels,       
                 sumDensityNodeLists = None,
                 useVelocityMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 evolveTotalEnergy = False,
                 gradhCorrection = False,
                 XSPH = True,
                 correctVelocityGradient = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 epsTensile = 0.0,
                 nTensile = 4.0,
                 damageRelieveRubble = True,
                 strengthInDamage = False,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()

        SolidFSISPHHydroBase%(dim)s.__init__(self,
                                          self._smoothingScaleMethod,
                                          dataBase,
                                          Q,
                                          slides,
                                          W,
                                          filter,
                                          cfl,
                                          surfaceForceCoefficient,
                                          densityStabilizationCoefficient, 
                                          specificThermalEnergyDiffusionCoefficient,
                                          xsphCoefficient,
                                          interfaceMethod,  
                                          kernelAveragingMethod,     
                                          sumDensityNodeLists,
                                          useVelocityMagnitudeForDt,
                                          compatibleEnergyEvolution,
                                          evolveTotalEnergy,
                                          gradhCorrection,
                                          XSPH,
                                          correctVelocityGradient,
                                          densityUpdate,
                                          HUpdate,
                                          epsTensile,
                                          nTensile,
                                          damageRelieveRubble,
                                          strengthInDamage,
                                          xmin,
                                          xmax)
        return
"""

for dim in dims:

    exec(SolidFSISPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidFSISPHHydro",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidFSISPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                       "classname"            : "SolidAFSISPHHydro",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})

def FSISPH(dataBase,
        W,
        Q = None,
        slides=None,
        filter = 0.0,
        cfl = 0.35,
        surfaceForceCoefficient=0.0,
        densityStabilizationCoefficient=0.1, 
        specificThermalEnergyDiffusionCoefficient=0.1, 
        xsphCoefficient=0.0,
        interfaceMethod=HLLCInterface, 
        kernelAveragingMethod = NeverAverageKernels,      
        sumDensityNodeLists=[],
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        correctVelocityGradient = True,    
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = True,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        RZ = False):

    # terms that are on deck or on their way out
    gradhCorrection = False
    densityUpdate = IntegrateDensity
    XSPH=False
    
    if compatibleEnergyEvolution and evolveTotalEnergy:
        raise RuntimeError, "compatibleEnergyEvolution and evolveTotalEnergy are incompatible"

    if gradhCorrection:
        raise RuntimeError, "gradhCorrection not implemented yet"

    if strengthInDamage and damageRelieveRubble:
        raise RuntimeError, "strengthInDamage and damageRelieveRubble are incompatible"

    # create the map nodelist --> index
    nodeLists = dataBase.nodeLists()
    nodeListMap = {}
    for i in range(dataBase.numNodeLists):          
        nodeListMap[nodeLists[i]]=i

    # for now use int 1/0 to indicate sum density or not
    sumDensitySwitch = [0]*dataBase.numNodeLists
    for nodeList in sumDensityNodeLists:
        i = nodeListMap[nodeList]
        sumDensitySwitch[i]=1
    sumDensitySwitch = vector_of_int(sumDensitySwitch)
            
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
        raise RuntimeError, "RZ is not implemented yet"
        # RZ ----------------------------------------
        '''
        if nsolid > 0:
            if ASPH:
                Constructor = SolidAFSISPHHydroRZ
            else:
                Constructor = SolidFSISPHHydroRZ
        else:
            if ASPH:
                Constructor = AFSISPHHydroRZ
            else:
                Constructor = FSISPHHydroRZ
        '''
    else:
    
        # Cartesian ---------------------------------
        if nsolid > 0:
            if ASPH:
                Constructor = eval("SolidAFSISPHHydro%id" % ndim)
            else:
                Constructor = eval("SolidFSISPHHydro%id" % ndim)
        else:
            raise RuntimeError, "currently only implemented for fluid nodelists"
            '''
            if ASPH:
                Constructor = eval("AFSISPHHydro%id" % ndim)
            else:
                Constructor = eval("FSISPHHydro%id" % ndim)
            '''

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 8.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("CRKSPHMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # slide surfaces.
    if not slides:
        contactTypes = vector_of_int([0]*(dataBase.numNodeLists**2))
        slides = eval("SlideSurface%id(contactTypes)" % (ndim))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "Q" : Q,
              "slides" : slides,
              "W" : W,
              "filter" : filter,
              "cfl" : cfl,
              "surfaceForceCoefficient" : surfaceForceCoefficient,
              "densityStabilizationCoefficient" : densityStabilizationCoefficient,      
              "specificThermalEnergyDiffusionCoefficient" : specificThermalEnergyDiffusionCoefficient,        
              "xsphCoefficient" : xsphCoefficient,
              "interfaceMethod" : interfaceMethod,
              "kernelAveragingMethod" : kernelAveragingMethod, 
              "sumDensityNodeLists" : sumDensitySwitch,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "gradhCorrection" : gradhCorrection,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    if nsolid > 0:
        kwargs.update({"damageRelieveRubble"      : damageRelieveRubble,
                       "strengthInDamage"         : strengthInDamage})

    # Build and return the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    result.slides = slides
    
    return result

#-------------------------------------------------------------------------------
# Provide a shorthand ASPH factory.
#-------------------------------------------------------------------------------
def AFSISPH(dataBase,
         W,
         Q = None,
         slides=None,
         filter = 0.0,
         cfl = 0.25,
         surfaceForceCoefficient = 0.0,
         densityStabilizationCoefficient = 0.0,  
         specificThermalEnergyDiffusionCoefficient = 0.0,
         xsphCoefficient = 0.0,   
         interfaceMethod = HLLCInterface,
         kernelAveragingMethod = NeverAverageKernels,                     
         sumDensityNodeLists = None,       
         useVelocityMagnitudeForDt = False,
         compatibleEnergyEvolution = True,
         evolveTotalEnergy = False,
         gradhCorrection = False,
         XSPH = False,
         correctVelocityGradient = True,
         densityUpdate = IntegrateDensity,
         HUpdate = IdealH,
         epsTensile = 0.0,
         nTensile = 4.0,
         damageRelieveRubble = False,
         strengthInDamage = False,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100)):
    return FSISPH(dataBase = dataBase,
               W = W,
               Q = Q,
               slides=slides,
               filter = filter,
               cfl = cfl,
               surfaceForceCoefficient = surfaceForceCoefficient,
               densityStabilizationCoefficient = densityStabilizationCoefficient, 
               specificThermalEnergyDiffusionCoefficient = specificThermalEnergyDiffusionCoefficient,  
               xsphCoefficient = xsphCoefficient,   
               interfaceMethod = interfaceMethod,    
               kernelAvragingMethod = kernelAveragingMethod,               
               sumDensityNodeLists = sumDensityNodeLists,          
               useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
               compatibleEnergyEvolution = compatibleEnergyEvolution,
               evolveTotalEnergy = evolveTotalEnergy,
               gradhCorrection = gradhCorrection,
               XSPH = XSPH,
               correctVelocityGradient = correctVelocityGradient,
               densityUpdate = densityUpdate,
               HUpdate = HUpdate,
               epsTensile = epsTensile,
               nTensile = nTensile,
               damageRelieveRubble = damageRelieveRubble,
               strengthInDamage = strengthInDamage,
               xmin = xmin,
               xmax = xmax,
               ASPH = True)
