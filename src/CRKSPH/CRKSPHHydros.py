from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(CRKSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 order,
                 filter,
                 cfl,
                 useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution,
                 evolveTotalEnergy,
                 XSPH,
                 densityUpdate,
                 HUpdate,
                 epsTensile,
                 nTensile):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        CRKSPHHydroBase%(dim)s.__init__(self,
                                        self._smoothingScaleMethod,
                                        dataBase,
                                        Q,
                                        order,
                                        filter,
                                        cfl,
                                        useVelocityMagnitudeForDt,
                                        compatibleEnergyEvolution,
                                        evolveTotalEnergy,
                                        XSPH,
                                        densityUpdate,
                                        HUpdate,
                                        epsTensile,
                                        nTensile)
        return
"""

#-------------------------------------------------------------------------------
# The generic SolidCRKSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidCRKSPHHydroFactoryString = """
class %(classname)s%(dim)s(SolidCRKSPHHydroBase%(dim)s):

    def __init__(self,
                 dataBase,
                 Q,
                 order,
                 filter,
                 cfl,
                 useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution,
                 evolveTotalEnergy,
                 XSPH,
                 densityUpdate,
                 HUpdate,
                 epsTensile,
                 nTensile,
                 damageRelieveRubble,
                 negativePressureInDamage):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        SolidCRKSPHHydroBase%(dim)s.__init__(self,
                                             self._smoothingScaleMethod,
                                             dataBase,
                                             Q,
                                             order,
                                             filter,
                                             cfl,
                                             useVelocityMagnitudeForDt,
                                             compatibleEnergyEvolution,
                                             evolveTotalEnergy,
                                             XSPH,
                                             densityUpdate,
                                             HUpdate,
                                             epsTensile,
                                             nTensile,
                                             damageRelieveRubble,
                                             negativePressureInDamage)
        return
"""

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
CRKSPHHydroRZFactoryString = """
class %(classname)s(CRKSPHHydroBaseRZ):

    def __init__(self,
                 dataBase,
                 Q,
                 order,
                 filter,
                 cfl,
                 useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution,
                 evolveTotalEnergy,
                 XSPH,
                 densityUpdate,
                 HUpdate,
                 epsTensile,
                 nTensile,
                 etaMinAxis):
        if GeometryRegistrar.coords() != CoordinateType.RZ:
            raise RuntimeError("Import from SpheralRZ before trying to use RZ physics")
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        CRKSPHHydroBaseRZ.__init__(self,
                                   self._smoothingScaleMethod,
                                   dataBase,
                                   Q,
                                   order,
                                   filter,
                                   cfl,
                                   useVelocityMagnitudeForDt,
                                   compatibleEnergyEvolution,
                                   evolveTotalEnergy,
                                   XSPH,
                                   densityUpdate,
                                   HUpdate,
                                   epsTensile,
                                   nTensile)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# The generic SolidCRKSPHHydro pattern.
#-------------------------------------------------------------------------------
SolidCRKSPHHydroRZFactoryString = """
class %(classname)s(SolidCRKSPHHydroBaseRZ):

    def __init__(self,
                 dataBase,
                 Q,
                 order,
                 filter,
                 cfl,
                 useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution,
                 evolveTotalEnergy,
                 XSPH,
                 densityUpdate,
                 HUpdate,
                 epsTensile,
                 nTensile,
                 damageRelieveRubble,
                 negativePressureInDamage,
                 etaMinAxis):
        if GeometryRegistrar.coords() != CoordinateType.RZ:
            raise RuntimeError("Import from SpheralRZ before trying to use RZ physics")
        self._smoothingScaleMethod = %(smoothingScaleMethod)s2d()
        SolidCRKSPHHydroBaseRZ.__init__(self,
                                        self._smoothingScaleMethod,
                                        dataBase,
                                        Q,
                                        order,
                                        filter,
                                        cfl,
                                        useVelocityMagnitudeForDt,
                                        compatibleEnergyEvolution,
                                        evolveTotalEnergy,
                                        XSPH,
                                        densityUpdate,
                                        HUpdate,
                                        epsTensile,
                                        nTensile,
                                        damageRelieveRubble,
                                        negativePressureInDamage)
        self.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        self.appendBoundary(self.zaxisBC)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "CRKSPHHydro",
                                     "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(CRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                     "classname"            : "ACRKSPHHydro",
                                     "smoothingScaleMethod" : "ASPHSmoothingScale"})

    exec(SolidCRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                          "classname"            : "SolidCRKSPHHydro",
                                          "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidCRKSPHHydroFactoryString % {"dim"                  : "%id" % dim,
                                          "classname"            : "SolidACRKSPHHydro",
                                          "smoothingScaleMethod" : "ASPHSmoothingScale"})

if 2 in dims:
    exec(CRKSPHHydroRZFactoryString % {"classname"            : "CRKSPHHydroRZ",
                                       "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(CRKSPHHydroRZFactoryString % {"classname"            : "ACRKSPHHydroRZ",
                                       "smoothingScaleMethod" : "ASPHSmoothingScale"})

    exec(SolidCRKSPHHydroRZFactoryString % {"classname"            : "SolidCRKSPHHydroRZ",
                                            "smoothingScaleMethod" : "SPHSmoothingScale"})
    exec(SolidCRKSPHHydroRZFactoryString % {"classname"            : "SolidACRKSPHHydroRZ",
                                            "smoothingScaleMethod" : "ASPHSmoothingScale"})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate CRKSPH hydro.
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
           negativePressureInDamage = False,
           ASPH = False,
           etaMinAxis = 0.1,
           crktype = "default"):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print "CRKSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported."
        print "             If you want some fluids active, provide SolidNodeList without a strength option specfied,"
        print "             which will result in fluid behaviour for those nodes."
        raise RuntimeError, "Cannot mix solid and fluid NodeLists."

    # Decide on the hydro object.
    if GeometryRegistrar.coords() == CoordinateType.RZ:

        # RZ ----------------------------------------
        if nsolid > 0:
            if ASPH:
                Constructor = SolidACRKSPHHydroRZ
            else:
                Constructor = SolidCRKSPHHydroRZ
        else:
            if ASPH:
                Constructor = ACRKSPHHydroRZ
            else:
                Constructor = CRKSPHHydroRZ

    else:

        crktype = crktype.lower()
        assert crktype in ("default", "variant")

        # Cartesian ---------------------------------
        if nsolid > 0:
            if ASPH:
                Constructor = eval("SolidACRKSPHHydro%id" % ndim)
            else:
                Constructor = eval("SolidCRKSPHHydro%id" % ndim)
        else:
            if ASPH:
                if crktype == "variant":
                    Constructor = eval("ACRKSPHVar%id" % ndim)
                else:
                    Constructor = eval("ACRKSPHHydro%id" % ndim)
            else:
                if crktype == "variant":
                    Constructor = eval("CRKSPHVar%id" % ndim)
                else:
                    Constructor = eval("CRKSPHHydro%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/4.0)
        Cq = 1.0*(dataBase.maxKernelExtent/4.0)**2
        Q = eval("CRKSPHMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Build the constructor arguments
    kwargs = {"dataBase" : dataBase,
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
        kwargs.update({"damageRelieveRubble" : damageRelieveRubble,
                       "negativePressureInDamage" : negativePressureInDamage})

    if GeometryRegistrar.coords() == CoordinateType.RZ:
        kwargs.update({"etaMinAxis" : etaMinAxis})

    # Build the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    return result

#-------------------------------------------------------------------------------
# ACRKSPH
#-------------------------------------------------------------------------------
def ACRKSPH(dataBase,
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
            negativePressureInDamage = False,
            etaMinAxis = 0.1):
    return CRKSPH(dataBase = dataBase,
                  Q = Q,
                  order = order,
                  filter = filter,
                  cfl = cfl,
                  useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                  compatibleEnergyEvolution = compatibleEnergyEvolution,
                  evolveTotalEnergy = evolveTotalEnergy,
                  XSPH = XSPH,
                  densityUpdate = densityUpdate,
                  HUpdate = HUpdate,
                  epsTensile = epsTensile,
                  nTensile = nTensile,
                  damageRelieveRubble = damageRelieveRubble,
                  negativePressureInDamage = negativePressureInDamage,
                  ASPH = True,
                  etaMinAxis = etaMinAxis)
