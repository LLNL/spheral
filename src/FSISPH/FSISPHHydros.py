from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

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
        damageRelieveRubble = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        RZ = False):
    ######################################################################
    # some of these parameters are inactive and possible on there was out.
    # strengthInDamage and damageRelieveRubble are old switches and are not
    # implemented in the code. RZ has not been implemented in FSISPH
    ######################################################################
    
    if compatibleEnergyEvolution and evolveTotalEnergy:
        raise RuntimeError("compatibleEnergyEvolution and evolveTotalEnergy are incompatible")

    if strengthInDamage and damageRelieveRubble:
        raise RuntimeError("strengthInDamage and damageRelieveRubble are incompatible")

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
        print("SPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("           If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("           which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Decide on the hydro object.
    if RZ:
        raise RuntimeError("RZ is not implemented yet")
    else:
        # Cartesian ---------------------------------
        if nsolid > 0:
            Constructor = eval("SolidFSISPHHydroBase%id" % ndim)
        else:
            raise RuntimeError("currently only implemented for solid nodelists")


    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 8.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # slide surfaces.
    if not slides:
        contactTypes = vector_of_int([0]*(dataBase.numNodeLists**2))
        slides = eval("SlideSurface%id(dataBase,contactTypes)" % ndim)

    # Smoothing scale update
    if ASPH:
        smoothingScaleMethod = eval("ASPHSmoothingScale%id()" % ndim)
    else:
        smoothingScaleMethod = eval("SPHSmoothingScale%id()" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"smoothingScaleMethod" : smoothingScaleMethod,
              "dataBase" : dataBase,
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
              "gradhCorrection" : False,
              "XSPH" : False,
              "correctVelocityGradient" : correctVelocityGradient,
              "densityUpdate" : IntegrateDensity,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "damageRelieveRubble" : damageRelieveRubble,
              "strengthInDamage" : strengthInDamage,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    result.slides = slides
    result._smoothingScaleMethod = smoothingScaleMethod
    
    return result

#-------------------------------------------------------------------------------
# Provide shorthand names for AFSISPH
#-------------------------------------------------------------------------------
def AFSISPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return SPH(*args, **kwargs)
