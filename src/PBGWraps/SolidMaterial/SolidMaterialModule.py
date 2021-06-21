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
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/SolidMaterialTypes.hh"' % srcdir)

        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        self.NinthOrderPolynomialFit = addObject(space, "NinthOrderPolynomialFit")

        for dim in self.dims:
            exec('''
EquationOfState%(dim)id = findObject(space, "EquationOfState%(dim)id")
Physics%(dim)id = findObject(space, "Physics%(dim)id")

self.SolidEquationOfState%(dim)id = addObject(space, "SolidEquationOfState%(dim)id", parent=EquationOfState%(dim)id, allow_subclassing=True)

self.PorousEquationOfState%(dim)id = addObject(space, "PorousEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.StrainPorosity%(dim)id = addObject(space, "StrainPorosity%(dim)id", parent=[Physics%(dim)id], allow_subclassing=True)

self.LinearPolynomialEquationOfState%(dim)id = addObject(space, "LinearPolynomialEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.GruneisenEquationOfState%(dim)id = addObject(space, "GruneisenEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.OsborneEquationOfState%(dim)id = addObject(space, "OsborneEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.TillotsonEquationOfState%(dim)id = addObject(space, "TillotsonEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)
self.MurnaghanEquationOfState%(dim)id = addObject(space, "MurnaghanEquationOfState%(dim)id", parent=self.SolidEquationOfState%(dim)id, allow_subclassing=True)

self.StrengthModel%(dim)id = addObject(space, "StrengthModel%(dim)id", allow_subclassing=True)
self.ConstantStrength%(dim)id = addObject(space, "ConstantStrength%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)
self.NullStrength%(dim)id = addObject(space, "NullStrength%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)
self.SteinbergGuinanStrength%(dim)id = addObject(space, "SteinbergGuinanStrength%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)
#self.SteinbergGuinanLundStrength%(dim)id = addObject(space, "SteinbergGuinanLundStrength%(dim)id", parent=self.SteinbergGuinanStrength%(dim)id, allow_subclassing=True)
self.JohnsonCookStrength%(dim)id = addObject(space, "JohnsonCookStrength%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)
self.CollinsStrength%(dim)id = addObject(space, "CollinsStrength%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)
self.PorousStrengthModel%(dim)id = addObject(space, "PorousStrengthModel%(dim)id", parent=self.StrengthModel%(dim)id, allow_subclassing=True)

self.PhysicsEvolvingMaterialLibrary%(dim)id = addObject(space, "PhysicsEvolvingMaterialLibrary%(dim)id", parent=[Physics%(dim)id, self.SolidEquationOfState%(dim)id, self.StrengthModel%(dim)id], allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        generateNinthOrderPolynomialFitBindings(self.NinthOrderPolynomialFit)

        for dim in self.dims:
            exec('''
generateSolidEquationOfStateBindings(self.SolidEquationOfState%(dim)id, %(dim)i)
generatePorousEquationOfStateBindings(self.PorousEquationOfState%(dim)id, %(dim)i)
generateStrainPorosityBindings(self.StrainPorosity%(dim)id, %(dim)i)

generateLinearPolynomialEquationOfStateBindings(self.LinearPolynomialEquationOfState%(dim)id, %(dim)i)
generateGruneisenEquationOfStateBindings(self.GruneisenEquationOfState%(dim)id, %(dim)i)
generateOsborneEquationOfStateBindings(self.OsborneEquationOfState%(dim)id, %(dim)i)
generateTillotsonEquationOfStateBindings(self.TillotsonEquationOfState%(dim)id, %(dim)i)
generateMurnaghanEquationOfStateBindings(self.MurnaghanEquationOfState%(dim)id, %(dim)i)

generateStrengthModelBindings(self.StrengthModel%(dim)id, %(dim)i)
generateConstantStrengthBindings(self.ConstantStrength%(dim)id, %(dim)i)
generateNullStrengthBindings(self.NullStrength%(dim)id, %(dim)i)
generateSteinbergGuinanStrengthBindings(self.SteinbergGuinanStrength%(dim)id, %(dim)i)
#generateSteinbergGuinanLundStrengthBindings(self.SteinbergGuinanLundStrength%(dim)id, %(dim)i)
generateJohnsonCookStrengthBindings(self.JohnsonCookStrength%(dim)id, %(dim)i)
generateCollinsStrengthBindings(self.CollinsStrength%(dim)id, %(dim)i)
generatePorousStrengthModelBindings(self.PorousStrengthModel%(dim)id, %(dim)i)

generatePhysicsEvolvingMaterialLibraryBindings(self.PhysicsEvolvingMaterialLibrary%(dim)id, %(dim)i)
''' % {"dim" : dim})

        return

#---------------------------------------------------------------------------
# SolidEquationOfState
#---------------------------------------------------------------------------
def generateSolidEquationOfStateBindings(x, ndim):

    # Constructors.
    x.add_constructor([param("double", "referenceDensity"),
                       param("double", "etamin"),
                       param("double", "etamax"),
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    # Methods.
    x.add_method("boundedEta", "double", [param("double", "rho")], is_const=True)
    x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

    # Attributes.
    x.add_instance_attribute("referenceDensity", "double", getter="referenceDensity", setter="referenceDensity")
    x.add_instance_attribute("etamin", "double", getter="etamin", setter="etamin")
    x.add_instance_attribute("etamax", "double", getter="etamax", setter="etamax")

    return

#---------------------------------------------------------------------------
# PorousEquationOfState
#---------------------------------------------------------------------------
def generatePorousEquationOfStateBindings(x, ndim):

    me = "Spheral::PorousEquationOfState%id" % ndim
    solideos = "Spheral::EquationOfState%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam(solideos, "solidEOS")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    const_ref_return_value(x, me, "%s::solidEOS" % me, solideos, [], "solidEOS")
    const_ref_return_value(x, me, "%s::alpha" % me, scalarfield, [], "alpha")

    # Attributes.
    x.add_instance_attribute("alpha0", "double", getter="alpha0", setter="alpha0")
    x.add_instance_attribute("c0", "double", getter="c0", setter="c0")

    return

#---------------------------------------------------------------------------
# LinearPolynomialEquationOfState
#---------------------------------------------------------------------------
def generateLinearPolynomialEquationOfStateBindings(x, ndim):

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
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                             param("double", "specificThermalEnergy")], is_const=True)
    x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
    x.add_method("pressure", "double", [param("double", "massDensity"),
                                        param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("temperature", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                     param("double", "temperature")],
                 is_const=True)
    x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                            param("double", "tempernature")],
                 is_const=True)
    x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                          param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("gamma", "double", [param("double", "massDensity"),
                                     param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("entropy", "double", [param("double", "massDensity"),
                                       param("double", "specificThermalEnergy")],
                 is_const=True)

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
def generateGruneisenEquationOfStateBindings(x, ndim):

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
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                             param("double", "specificThermalEnergy")], is_const=True)
    x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
    x.add_method("pressure", "double", [param("double", "massDensity"),
                                        param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("temperature", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                     param("double", "temperature")],
                 is_const=True)
    x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                            param("double", "tempernature")],
                 is_const=True)
    x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                          param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("gamma", "double", [param("double", "massDensity"),
                                     param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("entropy", "double", [param("double", "massDensity"),
                                       param("double", "specificThermalEnergy")],
                 is_const=True)

    # Attributes.
    x.add_instance_attribute("C0", "double", getter="C0", setter="C0")
    x.add_instance_attribute("S1", "double", getter="S1", setter="S1")
    x.add_instance_attribute("S2", "double", getter="S2", setter="S2")
    x.add_instance_attribute("S3", "double", getter="S3", setter="S3")
    x.add_instance_attribute("b", "double", getter="b", setter="b")
    x.add_instance_attribute("gamma0", "double", getter="gamma0", setter="gamma0")
    x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
    x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")
    x.add_instance_attribute("energyMultiplier", "double", getter="energyMultiplier", setter="energyMultiplier")

    return

#---------------------------------------------------------------------------
# OsborneEquationOfState
#---------------------------------------------------------------------------
def generateOsborneEquationOfStateBindings(x, ndim):

    # Constructors.
    x.add_constructor([param("double", "referenceDensity"),
                       param("double", "etamin"),
                       param("double", "etamax"),
                       param("double", "a1"),
                       param("double", "a2pos"),
                       param("double", "a2neg"),
                       param("double", "b0"),
                       param("double", "b1"),
                       param("double", "b2pos"),
                       param("double", "b2neg"),
                       param("double", "c0"),
                       param("double", "c1"),
                       param("double", "c2pos"),
                       param("double", "c2neg"),
                       param("double", "E0"),
                       param("double", "atomicWeight"),
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    x.add_method("DPDrho", "double", [param("double", "rho"),
                                      param("double", "specificThermalEnergy")], is_const=True)
    x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

    # Attributes.
    x.add_instance_attribute("a1", "double", getter="a1", setter="a1")
    x.add_instance_attribute("a2pos", "double", getter="a2pos", setter="a2pos")
    x.add_instance_attribute("a2neg", "double", getter="a2neg", setter="a2neg")
    x.add_instance_attribute("b0", "double", getter="b0", setter="b0")
    x.add_instance_attribute("b1", "double", getter="b1", setter="b1")
    x.add_instance_attribute("b2pos", "double", getter="b2pos", setter="b2pos")
    x.add_instance_attribute("b2neg", "double", getter="b2neg", setter="b2neg")
    x.add_instance_attribute("c0", "double", getter="c0", setter="c0")
    x.add_instance_attribute("c1", "double", getter="c1", setter="c1")
    x.add_instance_attribute("c2pos", "double", getter="c2pos", setter="c2pos")
    x.add_instance_attribute("c2neg", "double", getter="c2neg", setter="c2neg")
    x.add_instance_attribute("E0", "double", getter="E0", setter="E0")
    x.add_instance_attribute("atomicWeight", "double", getter="atomicWeight", setter="atomicWeight")
    x.add_instance_attribute("Cv", "double", getter="Cv", is_const=True)
    x.add_instance_attribute("externalPressure", "double", getter="externalPressure", setter="externalPressure")

    return

#---------------------------------------------------------------------------
# TillotsonEquationOfState
#---------------------------------------------------------------------------
def generateTillotsonEquationOfStateBindings(x, ndim):

    # Constructors.
    x.add_constructor([param("double", "referenceDensity"),
                       param("double", "etamin"),
                       param("double", "etamax"),
                       param("double", "etamin_solid"),
                       param("double", "etamax_solid"),
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
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

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

    x.add_method("pressure", "double", [param("double", "massDensity"),
                                        param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("temperature", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                     param("double", "temperature")],
                 is_const=True)
    x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                            param("double", "tempernature")],
                 is_const=True)
    x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                          param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("gamma", "double", [param("double", "massDensity"),
                                     param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("entropy", "double", [param("double", "massDensity"),
                                       param("double", "specificThermalEnergy")],
                 is_const=True)

    # Attributes.
    x.add_instance_attribute("etamin_solid", "double", getter="etamin_solid", setter="etamin_solid")
    x.add_instance_attribute("etamax_solid", "double", getter="etamax_solid", setter="etamax_solid")
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
# MurnaghanEquationOfState
#---------------------------------------------------------------------------
def generateMurnaghanEquationOfStateBindings(x, ndim):

    # Constructors.
    x.add_constructor([param("double", "referenceDensity"),
                       param("double", "etamin"),
                       param("double", "etamax"),
                       param("double", "n"),
                       param("double", "K"),
                       param("double", "atomicWeight"),
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "externalPressure", default_value="0.0"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    # Generic EOS interface.
    generateEquationOfStateVirtualBindings(x, ndim, False)

    # Methods.
    x.add_method("computeDPDrho", "double", [param("double", "rho"),
                                             param("double", "specificThermalEnergy")], is_const=True)
    x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
    x.add_method("pressure", "double", [param("double", "massDensity"),
                                        param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("temperature", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("specificThermalEnergy", "double", [param("double", "massDensity"),
                                                     param("double", "temperature")],
                 is_const=True)
    x.add_method("specificHeat", "double", [param("double", "massDensity"),
                                            param("double", "tempernature")],
                 is_const=True)
    x.add_method("soundSpeed", "double", [param("double", "massDensity"),
                                          param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("gamma", "double", [param("double", "massDensity"),
                                     param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("bulkModulus", "double", [param("double", "massDensity"),
                                           param("double", "specificThermalEnergy")],
                 is_const=True)
    x.add_method("entropy", "double", [param("double", "massDensity"),
                                       param("double", "specificThermalEnergy")],
                 is_const=True)

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
def generateStrengthModelVirtualBindings(x, ndim, pureVirtual):
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Methods.
    x.add_method("providesSoundSpeed", "bool", [], is_const=True, is_virtual=True)
    x.add_method("providesBulkModulus", "bool", [], is_const=True, is_virtual=True)
    x.add_method("shearModulus", None, [refparam(scalarfield, "shearModulus"),
                                        constrefparam(scalarfield, "density"),
                                        constrefparam(scalarfield, "specificThermalEnergy"),
                                        constrefparam(scalarfield, "pressure")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("yieldStrength", None, [refparam(scalarfield, "yieldStrength"),
                                         constrefparam(scalarfield, "density"),
                                         constrefparam(scalarfield, "specificThermalEnergy"),
                                         constrefparam(scalarfield, "pressure"),
                                         constrefparam(scalarfield, "plasticStrain"),
                                         constrefparam(scalarfield, "plasticStrainRate")],
                 is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("soundSpeed", None, [refparam(scalarfield, "soundSpeed"),
                                      constrefparam(scalarfield, "density"),
                                      constrefparam(scalarfield, "specificThermalEnergy"),
                                      constrefparam(scalarfield, "pressure"),
                                      constrefparam(scalarfield, "fluidSoundSpeed")],
                 is_const=True, is_virtual=True)
    x.add_method("bulkModulus", None, [refparam(scalarfield, "bulkModulus"),
                                       constrefparam(scalarfield, "density"),
                                       constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True)
    x.add_method("meltSpecificEnergy", None, [refparam(scalarfield, "meltSpecificEnergy"),
                                              constrefparam(scalarfield, "density"),
                                              constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True)
    x.add_method("coldSpecificEnergy", None, [refparam(scalarfield, "coldSpecificEnergy"),
                                              constrefparam(scalarfield, "density"),
                                              constrefparam(scalarfield, "specificThermalEnergy")],
                 is_const=True, is_virtual=True)
    return

#---------------------------------------------------------------------------
# StrengthModel
#---------------------------------------------------------------------------
def generateStrengthModelBindings(x, ndim):

    # Constructors.
    x.add_constructor([])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, True)

    return

#---------------------------------------------------------------------------
# ConstantStrength
#---------------------------------------------------------------------------
def generateConstantStrengthBindings(x, ndim):

    solidequationofstate = "Spheral::SolidEquationOfState%id" % ndim

    # Constructors.
    x.add_constructor([param("double", "mu0"),
                       param("double", "Y0")])
    x.add_constructor([param("double", "mu0"),
                       param("double", "Y0"),
                       constrefparam(solidequationofstate, "eos")])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Attributes.
    x.add_instance_attribute("mu0", "double", getter="mu0", is_const=True)
    x.add_instance_attribute("Y0", "double", getter="Y0", is_const=True)

    return

#---------------------------------------------------------------------------
# NullStrength
#---------------------------------------------------------------------------
def generateNullStrengthBindings(x, ndim):

    # Constructors.
    x.add_constructor([])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, False)

    return

#---------------------------------------------------------------------------
# NinthOrderPolynomialFit
#---------------------------------------------------------------------------
def generateNinthOrderPolynomialFitBindings(x):

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
    #x.add_constructor([constrefparam("Spheral::NinthOderPolynomialFit", "rhs")])

    # Methods.
    x.add_method("operator()", "double", [param("double", "x")], custom_name="__call__", is_const=True)

    return

#---------------------------------------------------------------------------
# SteinbergGuinanStrength
#---------------------------------------------------------------------------
def generateSteinbergGuinanStrengthBindings(x, ndim):

    solidequationofstate = "Spheral::SolidEquationOfState%id" % ndim
    ninthorderpolynomial = "Spheral::NinthOrderPolynomialFit"
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam(solidequationofstate, "eos"),
                       param("double", "G0"),
                       param("double", "Gmax"),
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
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Methods.
    x.add_method("meltAttenuation", "double", [param("double", "rho"), param("double", "eps")], is_const=True)
    x.add_method("computeTemperature", None, [refparam(scalarfield, "temperature"),
                                              constrefparam(scalarfield, "rho"), 
                                              constrefparam(scalarfield, "eps")], is_const=True)

    # Attributes.
    x.add_instance_attribute("G0", "double", getter="G0", is_const=True)
    x.add_instance_attribute("Gmax", "double", getter="G0", is_const=True)
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
# JohnsonCookStrength
#---------------------------------------------------------------------------
def generateJohnsonCookStrengthBindings(x, ndim):

    solidequationofstate = "Spheral::SolidEquationOfState%id" % ndim
    strengthmodel = "Spheral::StrengthModel%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam(solidequationofstate, "eos"),
                       constrefparam(strengthmodel, "shearModulusModel"),
                       param("double", "A"),
                       param("double", "B"),
                       param("double", "C"),
                       param("double", "C4"),
                       param("double", "m"),
                       param("double", "nhard"),
                       param("double", "epsdot0"),
                       param("double", "epsdotmin"),
                       param("double", "Tmelt"),
                       param("double", "Troom"),
                       param("double", "mu0", default_value="0.0"),
                       param("bool", "shearModulusScaling", default_value="false")])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Attributes.
    x.add_instance_attribute("A", "double", getter="A", is_const=True)
    x.add_instance_attribute("B", "double", getter="B", is_const=True)
    x.add_instance_attribute("C", "double", getter="C", is_const=True)
    x.add_instance_attribute("C4", "double", getter="C4", is_const=True)
    x.add_instance_attribute("m", "double", getter="m", is_const=True)
    x.add_instance_attribute("nhard", "double", getter="nhard", is_const=True)
    x.add_instance_attribute("epsdot0", "double", getter="epsdot0", is_const=True)
    x.add_instance_attribute("epsdotmin", "double", getter="epsdotmin", is_const=True)
    x.add_instance_attribute("Tmelt", "double", getter="Tmelt", is_const=True)
    x.add_instance_attribute("Troom", "double", getter="Troom", is_const=True)
    x.add_instance_attribute("mu0", "double", getter="mu0", is_const=True)
    x.add_instance_attribute("shearModulusScaling", "bool", getter="shearModulusScaling", is_const=True)

    return

#---------------------------------------------------------------------------
# CollinsStrength
#---------------------------------------------------------------------------
def generateCollinsStrengthBindings(x, ndim):

    solidequationofstate = "Spheral::SolidEquationOfState%id" % ndim
    strengthmodel = "Spheral::StrengthModel%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam(strengthmodel, "shearModulusModel"),
                       param("double", "mui"),
                       param("double", "Y0"),
                       param("double", "Ym")])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Attributes.
    x.add_instance_attribute("mui", "double", getter="mui", is_const=True)
    x.add_instance_attribute("Y0", "double", getter="Y0", is_const=True)
    x.add_instance_attribute("Ym", "double", getter="Ym", is_const=True)

    return

#---------------------------------------------------------------------------
# PorousStrengthModel
#---------------------------------------------------------------------------
def generatePorousStrengthModelBindings(x, ndim):

    me = "Spheral::PorousStrengthModel%id" % ndim
    strengthmodel = "Spheral::StrengthModel%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim

    # Constructors.
    x.add_constructor([constrefparam(strengthmodel, "solidStrength")])

    # Add the abstract interface.
    generateStrengthModelVirtualBindings(x, ndim, False)

    # Methods.
    const_ref_return_value(x, me, "%s::solidStrength" % me, strengthmodel, [], "solidStrength")
    const_ref_return_value(x, me, "%s::alpha" % me, scalarfield, [], "alpha")

    return

#---------------------------------------------------------------------------
# SteinbergGuinanLundStrength
#---------------------------------------------------------------------------
def generateSteinbergGuinanLundStrengthBindings(x, ndim):

    solidequationofstate = "Spheral::SolidEquationOfState%id" % ndim
    ninthorderpolynomial = "Spheral::NinthOrderPolynomialFit"

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
def generateStrainPorosityBindings(x, ndim):
    me = "Spheral::StrainPorosity%id" % ndim
    porouseos = "Spheral::PorousEquationOfState%id" % ndim
    porousstrength = "Spheral::PorousStrengthModel%id" % ndim
    nodelist = "Spheral::NodeList%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim
    vector_of_boundary = "vector_of_Boundary%id" % ndim
    fileio = "Spheral::FileIO"

    # Constructors.
    x.add_constructor([refparam(porouseos, "porousEOS"),
                       refparam(porousstrength, "porousStrength"),
                       constrefparam(nodelist, "nodeList"),
                       param("double", "phi0"),
                       param("double", "epsE"),
                       param("double", "epsX"),
                       param("double", "kappa"),
                       param("double", "gammaS0"),
                       param("double", "cS0"),
                       param("double", "c0")])

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
    x.add_instance_attribute("gammaS0", "double", getter="gammaS0", is_const=True)
    x.add_instance_attribute("cS0", "double", getter="cS0", is_const=True)
    x.add_instance_attribute("c0", "double", getter="c0", is_const=True)
    const_ref_return_value(x, me, "%s::porousEOS" % me, porouseos, [], "porousEOS")
    const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
    const_ref_return_value(x, me, "%s::alpha" % me, scalarfield, [], "alpha")
    const_ref_return_value(x, me, "%s::DalphaDt" % me, scalarfield, [], "DalphaDt")
    const_ref_return_value(x, me, "%s::strain" % me, scalarfield, [], "strain")
    const_ref_return_value(x, me, "%s::DstrainDt" % me, scalarfield, [], "DstrainDt")

    return

#---------------------------------------------------------------------------
# Geodyn
#---------------------------------------------------------------------------
def generatePhysicsEvolvingMaterialLibraryBindings(x, ndim):

    dim = "Spheral::Dim< %i >" % ndim
    me = "Spheral::PhysicsEvolvingMaterialLibrary%id" % ndim

    # Constructors.
    x.add_constructor([param("double", "referenceDensity"),
                       param("double", "etamin"),
                       param("double", "etamax"),
                       constrefparam("PhysicalConstants", "constants"),
                       param("double", "minimumPressure", default_value="-std::numeric_limits<double>::max()"),
                       param("double", "maximumPressure", default_value="std::numeric_limits<double>::max()"),
                       param("MaterialPressureMinType", "minPressureType", default_value="MaterialPressureMinType::PressureFloor")])

    return
