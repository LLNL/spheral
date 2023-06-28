from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
def SPH(W,
        WPi = None,
        WGrad = None,
        dataBase = None,
        Q = None,
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
        damageRelieveRubble = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        etaMinAxis = 0.1,
        ASPH = False,
        discretization = "Oslo"):

    # Check if we're running solid or fluid hydro
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("SPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("           If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("           which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    if GeometryRegistrar.coords() == CoordinateType.Spherical:
        assert ndim == 1
        if nsolid > 0:
            constructor = SolidSphericalSPHHydroOslo
        else:
            constructor = SphericalSPHHydroOslo

        # Build the spherical kernels
        print("Constructing Spherical kernels...")
        W3S1 = SphericalKernelOslo(W)
        W = W3S1
        if WPi:
            WPi3S1 = SphericalKernelOslo(WPi)
            WPi = WPi3S1
        if WGrad:
            WGrad3S1 = SphericalKernelOslo(WGrad)
            WGrad = WGrad3S1

    elif GeometryRegistrar.coords() == CoordinateType.RZ:
        assert ndim == 2
        if nsolid > 0:
            constructor = SolidSPHHydroBaseRZ
        else:
            constructor = SPHHydroBaseRZ

    else:
        if nsolid > 0:
            constructor = eval("SolidSPHHydroBase%id" % ndim)
        else:
            constructor = eval("SPHHydroBase%id" % ndim)

    # Fill out the set of kernels
    if WPi is None:
        WPi = W
    if WGrad is None:
        WGrad = W

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

    # Smoothing scale update
    if ASPH:
        smoothingScaleMethod = eval("ASPHSmoothingScale%id()" % ndim)
    else:
        smoothingScaleMethod = eval("SPHSmoothingScale%id()" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"smoothingScaleMethod" : smoothingScaleMethod,
              "W" : W,
              "WPi" : WPi,
              "dataBase" : dataBase,
              "Q" : Q,
              "filter" : filter,
              "cfl" : cfl,
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
        kwargs.update({"WGrad"                    : WGrad,
                       "damageRelieveRubble"      : damageRelieveRubble,
                       "strengthInDamage"         : strengthInDamage})

    # Build the SPH hydro
    result = constructor(**kwargs)
    result.Q = Q
    result._smoothingScaleMethod = smoothingScaleMethod

    # In spherical coordinates, preserve our locally constructed spherical kernels
    # and add the origin enforcement boundary
    if GeometryRegistrar.coords() == CoordinateType.Spherical:
        result.originBC = SphericalOriginBoundary()
        result.appendBoundary(result.originBC)
        result._W3S1 = W
        result._WPi3S1 = WPi
        result._WGrad3S1 = WGrad

    # If we're using area-weighted RZ, we need to reflect from the axis
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        result.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        result.appendBoundary(result.zaxisBC)

    return result

#-------------------------------------------------------------------------------
# Provide shorthand names for ASPH
#-------------------------------------------------------------------------------
def ASPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return SPH(*args, **kwargs)
