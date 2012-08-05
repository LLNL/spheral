from pybindgen import *

from ref_return_value import *
from MaterialModule import  generateEquationOfStateVirtualBindings
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SolidMaterial:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"SolidMaterial/SolidMaterialTypes.hh"')

        # Namespace.
        SolidSpheral = mod.add_cpp_namespace("Spheral")
        space = SolidSpheral.add_cpp_namespace("SolidMaterial")

        Spheral = mod.add_cpp_namespace("Spheral")
        Material = Spheral.add_cpp_namespace("Material")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")

        self.unitSet = ("CGS", "MKS", "Solar")
        self.dimSet = (1, 2, 3)

        self.StrengthModel = addObject(space, "StrengthModel", allow_subclassing=True)
        self.ConstantStrength = addObject(space, "ConstantStrength", parent=self.StrengthModel, allow_subclassing=True)
        self.NullStrength = addObject(space, "NullStrength", parent=self.StrengthModel, allow_subclassing=True)
        self.NinthOrderPolynomialFit = addObject(space, "NinthOrderPolynomialFit")

        for dim in self.dimSet:
            exec('''
EquationOfState%(dim)id = Material.wrapObjs["EquationOfState%(dim)id"]
Physics%(dim)id = PhysicsSpace.wrapObjs["Physics%(dim)id"]
self.SolidEquationOfState%(dim)id = addObject(space, "SolidEquationOfState%(dim)id", parent=EquationOfState%(dim)id, allow_subclassing=True)
self.PorousEquationOfState%(dim)id = addObject(space, "PorousEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.StrainPorosity%(dim)id = addObject(space, "StrainPorosity%(dim)id", parent=[Physics%(dim)id], allow_subclassing=True)
''' % {"dim" : dim})
            for units in self.unitSet:
                exec('''
self.LinearPolynomialEquationOfState%(units)s%(dim)id = addObject(space, "LinearPolynomialEquationOfState%(units)s%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.GruneisenEquationOfState%(units)s%(dim)id = addObject(space, "GruneisenEquationOfState%(units)s%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.TillotsonEquationOfState%(units)s%(dim)id = addObject(space, "TillotsonEquationOfState%(units)s%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.MurnahanEquationOfState%(units)s%(dim)id = addObject(space, "MurnahanEquationOfState%(units)s%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.SteinbergGuinanStrength%(units)s%(dim)id = addObject(space, "SteinbergGuinanStrength%(units)s%(dim)id", parent=self.StrengthModel, allow_subclassing=True)
self.SteinbergGuinanLundStrength%(units)s%(dim)id = addObject(space, "SteinbergGuinanLundStrength%(units)s%(dim)id", parent=self.SteinbergGuinanStrength%(units)s%(dim)id, allow_subclassing=True)
''' % {"dim" : dim,
       "units" : units})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateStrengthModelBindings(self.StrengthModel)
        self.generateConstantStrengthBindings(self.ConstantStrength)
        self.generateNullStrengthBindings(self.NullStrength)
        self.generateNinthOrderPolynomialFitBindings(self.NinthOrderPolynomialFit)

        for dim in self.dimSet:
            exec('''
self.generateSolidEquationOfStateBindings(self.SolidEquationOfState%(dim)id, %(dim)i)
self.generatePorousEquationOfStateBindings(self.PorousEquationOfState%(dim)id, %(dim)i)
self.generateStrainPorosityBindings(self.StrainPorosity%(dim)id, %(dim)i)
''' % {"dim" : dim})

            for units in self.unitSet:
                exec('''
self.generateLinearPolynomialEquationOfStateBindings(self.LinearPolynomialEquationOfState%(units)s%(dim)id, %(dim)i)
self.generateGruneisenEquationOfStateBindings(self.GruneisenEquationOfState%(units)s%(dim)id, %(dim)i)
self.generateTillotsonEquationOfStateBindings(self.TillotsonEquationOfState%(units)s%(dim)id, %(dim)i)
self.generateMurnahanEquationOfStateBindings(self.MurnahanEquationOfState%(units)s%(dim)id, %(dim)i)

self.generateSteinbergGuinanStrengthBindings(self.SteinbergGuinanStrength%(units)s%(dim)id, %(dim)i)
self.generateSteinbergGuinanLundStrengthBindings(self.SteinbergGuinanLundStrength%(units)s%(dim)id, %(dim)i)
''' % {"dim" : dim,
       "units" : units})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # SolidEquationOfState
    #---------------------------------------------------------------------------
    def generateSolidEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "referenceDensity"),
                           param("double", "etamin"),
                           param("double", "etamax"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Methods.
        x.add_method("boundedEta", "double", [param("double", "rho")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("referenceDensity", "double", getter="referenceDensity", setter="referenceDensity")
        x.add_instance_attribute("etamin", "double", getter="etamin", setter="etamin")
        x.add_instance_attribute("etamax", "double", getter="etamax", setter="etamax")

        return

    #---------------------------------------------------------------------------
    # SolidEquationOfState
    #---------------------------------------------------------------------------
    def generatePorousEquationOfStateBindings(self, x, ndim):

        me = "Spheral::SolidMaterial::PorousEquationOfState%id" % ndim
        solideos = "Spheral::SolidMaterial::SolidEquationOfState%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(solideos, "solidEOS")])

        # Generic EOS interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        # Methods.
        const_ref_return_value(x, me, "%s::solidEOS" % me, solideos, [], "solidEOS")
        const_ref_return_value(x, me, "%s::alpha" % me, scalarfield, [], "alpha")

        return

    #---------------------------------------------------------------------------
    # LinearPolynomialEquationOfState
    #---------------------------------------------------------------------------
    def generateLinearPolynomialEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "referenceDensity"),
                           param("double", "etamin"),
                           param("double", "etamax"),
                           param("double", "a0"),
                           param("double", "a1"),
                           param("double", "a2"),
                           param("double", "a3"),
                           param("double", "b0"),
                           param("double", "b1"),
                           param("double", "b2"),
                           param("double", "atomicWeight"),
                           param("double", "externalPressure", default_value="0.0"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Generic EOS interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                                 param("double", "specificThermalEnergy")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("a0", "double", getter="a0", setter="a0")
        x.add_instance_attribute("a1", "double", getter="a1", setter="a1")
        x.add_instance_attribute("a2", "double", getter="a2", setter="a2")
        x.add_instance_attribute("a3", "double", getter="a3", setter="a3")
        x.add_instance_attribute("b0", "double", getter="b0", setter="b0")
        x.add_instance_attribute("b1", "double", getter="b1", setter="b1")
        x.add_instance_attribute("b2", "double", getter="b2", setter="b2")
        x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")
        x.add_instance_attribute("minimumPressure", "double", getter="minimumPressure", setter="minimumPressure")

        return

    #---------------------------------------------------------------------------
    # GruneisenEquationOfState
    #---------------------------------------------------------------------------
    def generateGruneisenEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "referenceDensity"),
                           param("double", "etamin"),
                           param("double", "etamax"),
                           param("double", "C0"),
                           param("double", "S1"),
                           param("double", "S2"),
                           param("double", "S3"),
                           param("double", "gamma0"),
                           param("double", "b"),
                           param("double", "atomicWeight"),
                           param("double", "externalPressure", default_value="0.0"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Generic EOS interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                                 param("double", "specificThermalEnergy")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("C0", "double", getter="C0", setter="C0")
        x.add_instance_attribute("S1", "double", getter="S1", setter="S1")
        x.add_instance_attribute("S2", "double", getter="S2", setter="S2")
        x.add_instance_attribute("S3", "double", getter="S3", setter="S3")
        x.add_instance_attribute("b", "double", getter="b", setter="b")
        x.add_instance_attribute("gamma0", "double", getter="gamma0", setter="gamma0")
        x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")

        return

    #---------------------------------------------------------------------------
    # TillotsonEquationOfState
    #---------------------------------------------------------------------------
    def generateTillotsonEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "referenceDensity"),
                           param("double", "etamin"),
                           param("double", "etamax"),
                           param("double", "a"),
                           param("double", "b"),
                           param("double", "A"),
                           param("double", "B"),
                           param("double", "alpha"),
                           param("double", "beta"),
                           param("double", "eps0"),
                           param("double", "epsLiquid"),
                           param("double", "epsVapor"),
                           param("double", "atomicWeight"),
                           param("double", "externalPressure", default_value="0.0"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Generic EOS interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                                 param("double", "specificThermalEnergy")], is_const=True)
        x.add_method("computePhi", "double", [param("double", "eta"), param("double", "eps")], is_const=True)
        x.add_method("computeP1", "double", [param("double", "mu"), param("double", "P2")], is_const=True)
        x.add_method("computeP2", "double", [param("double", "phi"), param("double", "mu"), param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("computeP4", "double", [param("double", "phi"), param("double", "mu"), param("double", "eta"), param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("compute_dphidrho_eps", "double", [param("double", "rho0"), param("double", "eta"), param("double", "eps")], is_const=True)
        x.add_method("compute_dP1drho_eps", "double", [param("double", "rho0"), param("double", "mu"), param("double", "dP2drho_eps")], is_const=True)
        x.add_method("compute_dP2drho_eps", "double", [param("double", "phi"), param("double", "dphidrho_eps"), param("double", "rho0"), param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("compute_dP4drho_eps", "double", [param("double", "phi"), param("double", "dphidrho_eps"), param("double", "rho0"), param("double", "eta"), param("double", "mu"), param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("compute_dphideps_rho", "double", [param("double", "eta"), param("double", "eps")], is_const=True)
        x.add_method("compute_dP2deps_rho", "double", [param("double", "phi"), param("double", "dphideps_rho"), param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("compute_dP4deps_rho", "double", [param("double", "phi"), param("double", "dphideps_rho"), param("double", "eta"), param("double", "rho"), param("double", "eps")], is_const=True)

        # Attributes.
        x.add_instance_attribute("a", "double", getter="a", setter="a")
        x.add_instance_attribute("b", "double", getter="b", setter="b")
        x.add_instance_attribute("A", "double", getter="A", setter="A")
        x.add_instance_attribute("B", "double", getter="B", setter="B")
        x.add_instance_attribute("alpha", "double", getter="alpha", setter="alpha")
        x.add_instance_attribute("beta", "double", getter="beta", setter="beta")
        x.add_instance_attribute("eps0", "double", getter="eps0", setter="eps0")
        x.add_instance_attribute("epsLiquid", "double", getter="epsLiquid", setter="epsLiquid")
        x.add_instance_attribute("epsVapor", "double", getter="epsVapor", setter="epsVapor")
        x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")

        return

    #---------------------------------------------------------------------------
    # MurnahanEquationOfState
    #---------------------------------------------------------------------------
    def generateMurnahanEquationOfStateBindings(self, x, ndim):

        # Constructors.
        x.add_constructor([param("double", "referenceDensity"),
                           param("double", "etamin"),
                           param("double", "etamax"),
                           param("double", "n"),
                           param("double", "K"),
                           param("double", "atomicWeight"),
                           param("double", "externalPressure", default_value="0.0"),
                           param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                           param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()")])

        # Generic EOS interface.
        generateEquationOfStateVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                                 param("double", "specificThermalEnergy")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("n", "double", getter="n", setter="n")
        x.add_instance_attribute("K", "double", getter="K", setter="K")
        x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
        x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")
        x.add_instance_attribute("minimumPressure", "double", getter="minimumPressure", setter="minimumPressure")

        return

    #---------------------------------------------------------------------------
    # StrengthModel (virtual interface)
    #---------------------------------------------------------------------------
    def generateStrengthModelVirtualBindings(self, x, pureVirtual):

        # Methods.
        x.add_method("shearModulus", "double", [param("double", "density"),
                                                param("double", "specificThermalEnergy"),
                                                param("double", "pressure")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("yieldStrength", "double", [param("double", "density"),
                                                 param("double", "specificThermalEnergy"),
                                                 param("double", "pressure"),
                                                 param("double", "plasticStrain"),
                                                 param("double", "plasticStrainRate")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("soundSpeed", "double", [param("double", "density"),
                                              param("double", "specificThermalEnergy"),
                                              param("double", "pressure"),
                                              param("double", "fluidSoundSpeed")],
                     is_const=True, is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # StrengthModel
    #---------------------------------------------------------------------------
    def generateStrengthModelBindings(self, x):

        # Constructors.
        x.add_constructor([])

        # Add the abstract interface.
        self.generateStrengthModelVirtualBindings(x, True)

        return

    #---------------------------------------------------------------------------
    # ConstantStrength
    #---------------------------------------------------------------------------
    def generateConstantStrengthBindings(self, x):

        # Constructors.
        x.add_constructor([param("double", "mu0"),
                           param("double", "Y0")])

        # Add the abstract interface.
        self.generateStrengthModelVirtualBindings(x, False)

        return

    #---------------------------------------------------------------------------
    # NullStrength
    #---------------------------------------------------------------------------
    def generateNullStrengthBindings(self, x):

        # Constructors.
        x.add_constructor([])

        # Add the abstract interface.
        self.generateStrengthModelVirtualBindings(x, False)

        return

    #---------------------------------------------------------------------------
    # NinthOrderPolynomialFit
    #---------------------------------------------------------------------------
    def generateNinthOrderPolynomialFitBindings(self, x):

        # Constructors.
        x.add_constructor([param("double", "C0"),
                           param("double", "C1"),
                           param("double", "C2"),
                           param("double", "C3"),
                           param("double", "C4"),
                           param("double", "C5"),
                           param("double", "C6"),
                           param("double", "C7"),
                           param("double", "C8"),
                           param("double", "C9")])
        #x.add_constructor([constrefparam("Spheral::SolidMaterial::NinthOderPolynomialFit", "rhs")])

        # Methods.
        x.add_method("operator()", "double", [param("double", "x")], custom_name="__call__", is_const=True)

        return

    #---------------------------------------------------------------------------
    # SteinbergGuinanStrength
    #---------------------------------------------------------------------------
    def generateSteinbergGuinanStrengthBindings(self, x, ndim):

        solidequationofstate = "Spheral::SolidMaterial::SolidEquationOfState%id" % ndim
        ninthorderpolynomial = "Spheral::SolidMaterial::NinthOrderPolynomialFit"

        # Constructors.
        x.add_constructor([constrefparam(solidequationofstate, "eos"),
                           param("double", "G0"),
                           param("double", "A"),
                           param("double", "B"),
                           param("double", "Y0"),
                           param("double", "Ymax"),
                           param("double", "Yp"),
                           param("double", "beta"),
                           param("double", "gamma0"),
                           param("double", "nhard"),
                           constrefparam(ninthorderpolynomial, "coldEnergyFit"),
                           constrefparam(ninthorderpolynomial, "meltEnergyFit")])

        # Add the abstract interface.
        self.generateStrengthModelVirtualBindings(x, False)

        # Methods.
        x.add_method("meltAttenuation", "double", [param("double", "rho"), param("double", "eps")], is_const=True)
        x.add_method("computeTemperature", "double", [param("double", "rho"), param("double", "eps")], is_const=True)

        # Attributes.
        x.add_instance_attribute("G0", "double", getter="G0", is_const=True)
        x.add_instance_attribute("A", "double", getter="A", is_const=True)
        x.add_instance_attribute("B", "double", getter="B", is_const=True)
        x.add_instance_attribute("Y0", "double", getter="Y0", is_const=True)
        x.add_instance_attribute("Ymax", "double", getter="Ymax", is_const=True)
        x.add_instance_attribute("Yp", "double", getter="Yp", is_const=True)
        x.add_instance_attribute("beta", "double", getter="beta", is_const=True)
        x.add_instance_attribute("gamma0", "double", getter="gamma0", is_const=True)
        x.add_instance_attribute("nhard", "double", getter="nhard", is_const=True)
        x.add_instance_attribute("coldEnergyFit", ninthorderpolynomial, getter="coldEnergyFit", is_const=True)
        x.add_instance_attribute("meltEnergyFit", ninthorderpolynomial, getter="meltEnergyFit", is_const=True)

        return

    #---------------------------------------------------------------------------
    # SteinbergGuinanLundStrength
    #---------------------------------------------------------------------------
    def generateSteinbergGuinanLundStrengthBindings(self, x, ndim):

        solidequationofstate = "Spheral::SolidMaterial::SolidEquationOfState%id" % ndim
        ninthorderpolynomial = "Spheral::SolidMaterial::NinthOrderPolynomialFit"

        # Constructors.
        x.add_constructor([constrefparam(solidequationofstate, "eos"),
                           param("double", "G0"),
                           param("double", "A"),
                           param("double", "B"),
                           param("double", "Y0"),
                           param("double", "Ymax"),
                           param("double", "Yp"),
                           param("double", "beta"),
                           param("double", "gamma0"),
                           param("double", "nhard"),
                           param("double", "C1"),
                           param("double", "C2"),
                           param("double", "UK"),
                           param("double", "YP"),
                           param("double", "YTmax"),
                           constrefparam(ninthorderpolynomial, "coldEnergyFit"),
                           constrefparam(ninthorderpolynomial, "meltEnergyFit")])

        # Methods.
        x.add_method("yieldStrength", "double", [param("double", "density"),
                                                 param("double", "specificThermalEnergy"),
                                                 param("double", "pressure"),
                                                 param("double", "plasticStrain"),
                                                 param("double", "plasticStrainRate")],
                     is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("C1", "double", getter="C1", is_const=True)
        x.add_instance_attribute("C2", "double", getter="C2", is_const=True)
        x.add_instance_attribute("UK", "double", getter="UK", is_const=True)
        x.add_instance_attribute("YP", "double", getter="YP", is_const=True)
        x.add_instance_attribute("YTmax", "double", getter="YTmax", is_const=True)

        return

    #---------------------------------------------------------------------------
    # StrainPorosity
    #---------------------------------------------------------------------------
    def generateStrainPorosityBindings(self, x, ndim):
        me = "Spheral::SolidMaterial::StrainPorosity%id" % ndim
        porouseos = "Spheral::SolidMaterial::PorousEquationOfState%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([refparam(porouseos, "porousEOS"),
                           constrefparam(nodelist, "nodeList"),
                           param("double", "phi0"),
                           param("double", "epsE"),
                           param("double", "epsX"),
                           param("double", "kappa")])

        # Physics interface.
        generatePhysicsVirtualBindings(x, ndim, False)

        # Methods.
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, 
                     [refparam(fileio, "file"), refparam("std::string", "pathName")],
                     is_const=True, is_virtual=True)
        x.add_method("restoreState", None, 
                     [constrefparam(fileio, "file"), refparam("std::string", "pathName")],
                     is_virtual=True)
        
        # Attributes.
        x.add_instance_attribute("phi0", "double", getter="phi0", is_const=True)
        x.add_instance_attribute("alpha0", "double", getter="alpha0", is_const=True)
        x.add_instance_attribute("epsE", "double", getter="epsE", is_const=True)
        x.add_instance_attribute("epsX", "double", getter="epsX", is_const=True)
        x.add_instance_attribute("epsC", "double", getter="epsC", is_const=True)
        x.add_instance_attribute("kappa", "double", getter="kappa", is_const=True)
        const_ref_return_value(x, me, "%s::porousEOS" % me, porouseos, [], "porousEOS")
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
        const_ref_return_value(x, me, "%s::alpha" % me, scalarfield, [], "alpha")
        const_ref_return_value(x, me, "%s::DalphaDt" % me, scalarfield, [], "DalphaDt")
        const_ref_return_value(x, me, "%s::strain" % me, scalarfield, [], "strain")
        const_ref_return_value(x, me, "%s::DstrainDt" % me, scalarfield, [], "DstrainDt")

        return

