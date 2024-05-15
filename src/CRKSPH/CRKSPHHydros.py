from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
def CRKSPH(dataBase,
           Q = None,
           order = RKOrder.LinearOrder,
           filter = 0.0,
           cfl = 0.25,
           useVelocityMagnitudeForDt = False,
           compatibleEnergyEvolution = True,
           evolveTotalEnergy = False,
           XSPH = True,
           densityUpdate = RigorousSumDensity,
           HUpdate = IdealH,
           epsTensile = 0.0,
           nTensile = 4.0,
           damageRelieveRubble = False,
           ASPH = False,
           etaMinAxis = 0.1,
           crktype = "default",
           smoothingScaleMethod = None):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("CRKSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("             If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("             which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        # RZ ----------------------------------------
        assert ndim == 2
        if nsolid > 0:
            constructor = SolidCRKSPHHydroBaseRZ
        else:
            constructor = CRKSPHHydroBaseRZ
    else:
        # Cartesian ---------------------------------
        crktype = crktype.lower()
        assert crktype in ("default", "variant")
        if nsolid > 0:
            constructor = eval("SolidCRKSPHHydroBase%id" % ndim)
        else:
            if crktype == "variant":
                constructor = eval("CRKSPHVariant%id" % ndim)
            else:
                constructor = eval("CRKSPHHydroBase%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/4.0)
        Cq = 1.0*(dataBase.maxKernelExtent/4.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Smoothing scale update
    if smoothingScaleMethod is None:
        if ASPH:
            smoothingScaleMethod = eval("ASPHSmoothingScale%id()" % ndim)
        else:
            smoothingScaleMethod = eval("SPHSmoothingScale%id()" % ndim)

    # Build the constructor arguments
    kwargs = {"smoothingScaleMethod" : smoothingScaleMethod,
              "dataBase" : dataBase,
              "Q" : Q,
              "order" : order,
              "filter" : filter,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "densityUpdate" : densityUpdate,
              "HUpdate" : HUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile}

    if nsolid > 0:
        kwargs.update({"damageRelieveRubble" : damageRelieveRubble})

    if GeometryRegistrar.coords() == CoordinateType.RZ:
        kwargs.update({"etaMinAxis" : etaMinAxis})

    # Build the thing.
    result = constructor(**kwargs)
    result.Q = Q
    result._smoothingScaleMethod = smoothingScaleMethod
    return result

#-------------------------------------------------------------------------------
# ACRKSPH
#-------------------------------------------------------------------------------
def ACRKSPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return CRKSPH(*args, **kwargs)
