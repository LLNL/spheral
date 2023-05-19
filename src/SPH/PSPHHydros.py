from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic PSPHHydro pattern.
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
         HopkinsConductivity = False,
         sumMassDensityOverAllNodeLists = True,
         densityUpdate = RigorousSumDensity,
         HUpdate = IdealH,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100),
         ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("PSPH Warning: you have provided solid NodeLists, but PSPH currently does not have a solid option.")
        print("              The fluid limit will be provided for now.")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    constructor = eval("PSPHHydroBase%id" % ndim)

    # Fill out the set of kernels
    if WPi is None:
        WPi = W

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 2.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

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
              "W" : W,
              "WPi" : WPi,
              "filter" : filter,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "HopkinsConductivity" : HopkinsConductivity,
              "sumMassDensityOverAllNodeLists" : sumMassDensityOverAllNodeLists,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build the thing
    result = constructor(**kwargs)
    result.Q = Q
    result._smoothingScaleMethod = smoothingScaleMethod
    
    return result

#-------------------------------------------------------------------------------
# Provide a shorthand APSPH factory.
#-------------------------------------------------------------------------------
def PASPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return PSPH(*args, **kwargs)
