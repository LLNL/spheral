from pybindgen import *

from ref_return_value import *
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SolidSPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/SolidSPHTypes.hh"' % srcdir)
    
        # Namespaces.
        SolidSpheral = mod.add_cpp_namespace("Spheral")
        space = SolidSpheral.add_cpp_namespace("SolidSPHSpace")

        Spheral = mod.add_cpp_namespace("Spheral")
        SPHSpace = Spheral.add_cpp_namespace("SPHSpace")

        sphhydrobase1d = findObject(SPHSpace, "SPHHydroBase1d")
        sphhydrobase2d = findObject(SPHSpace, "SPHHydroBase2d")
        sphhydrobase3d = findObject(SPHSpace, "SPHHydroBase3d")

        # Expose types.
        self.NodeCoupling = addObject(Spheral, "NodeCoupling", allow_subclassing=True)

        self.DamagedNodeCoupling1d = addObject(Spheral, "DamagedNodeCoupling1d", allow_subclassing=True, parent=self.NodeCoupling)
        self.DamagedNodeCoupling2d = addObject(Spheral, "DamagedNodeCoupling2d", allow_subclassing=True, parent=self.NodeCoupling)
        self.DamagedNodeCoupling3d = addObject(Spheral, "DamagedNodeCoupling3d", allow_subclassing=True, parent=self.NodeCoupling)

        self.SolidSPHHydroBase1d = addObject(space, "SolidSPHHydroBase1d", allow_subclassing=True, parent=sphhydrobase1d)
        self.SolidSPHHydroBase2d = addObject(space, "SolidSPHHydroBase2d", allow_subclassing=True, parent=sphhydrobase2d)
        self.SolidSPHHydroBase3d = addObject(space, "SolidSPHHydroBase3d", allow_subclassing=True, parent=sphhydrobase3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateNodeCouplingBindings(self.NodeCoupling)

        self.generateDamagedNodeCouplingBindings(self.DamagedNodeCoupling1d, 1)
        self.generateDamagedNodeCouplingBindings(self.DamagedNodeCoupling2d, 2)
        self.generateDamagedNodeCouplingBindings(self.DamagedNodeCoupling3d, 3)

        self.generateSolidSPHHydroBaseBindings(self.SolidSPHHydroBase1d, 1)
        self.generateSolidSPHHydroBaseBindings(self.SolidSPHHydroBase2d, 2)
        self.generateSolidSPHHydroBaseBindings(self.SolidSPHHydroBase3d, 3)
        
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["SolidSPHSpace"]

    #---------------------------------------------------------------------------
    # Bindings (NodeCoupling).
    #---------------------------------------------------------------------------
    def generateNodeCouplingBindings(self, x):

        # Object names.
        me = "Spheral::NodeCoupling"

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("operator()", "double", [param("unsigned int", "nodeListi"),
                                              param("unsigned int", "i"),
                                              param("unsigned int", "nodeListj"),
                                              param("unsigned int", "j")],
                     custom_name = "__call__",
                     is_const = True,
                     is_virtual = True)
        return

    #---------------------------------------------------------------------------
    # Bindings (DamagedNodeCoupling).
    #---------------------------------------------------------------------------
    def generateDamagedNodeCouplingBindings(self, x, ndim):

        # Object names.
        me = "Spheral::DamagedNodeCoupling%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(symtensorfieldlist, "damage"),
                           constrefparam(vectorfieldlist, "damageGradient"),
                           constrefparam(symtensorfieldlist, "H")])

        return

    #---------------------------------------------------------------------------
    # Bindings (SolidSPHHydroBase).
    #---------------------------------------------------------------------------
    def generateSolidSPHHydroBaseBindings(self, x, ndim):

        # Object names.
        me = "Spheral::SolidSPHSpace::SolidSPHHydroBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                           constrefparam(tablekernel, "W"),
                           constrefparam(tablekernel, "WPi"),
                           refparam(artificialviscosity, "Q"),
                           param("double", "filter", default_value="0.0"),
                           param("double", "cfl", default_value="0.5"),
                           param("int", "useVelocityMagnitudeForDt", default_value="false"),
                           param("int", "compatibleEnergyEvolution", default_value="true"),
                           param("int", "gradhCorrection", default_value="false"),
                           param("int", "XSPH", default_value="true"),
                           param("int", "correctVelocityGradient", default_value="true"),
                           param("int", "sumMassDensityOverAllNodeLists", default_value="false"),
                           param("MassDensityType", "densityUpdate", default_value="Spheral::PhysicsSpace::RigorousSumDensity"),
                           param("HEvolutionType", "HUpdate", default_value="Spheral::PhysicsSpace::IdealH"),
                           param("double", "epsTensile", default_value="0.3"),
                           param("double", "nTensile", default_value="4.0"),
                           param(vector, "xmin", default_value="%s(-1e10, -1e10, -1e10)" % vector),
                           param(vector, "xmax", default_value="%s( 1e10,  1e10,  1e10)" % vector)])

        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        x.add_method("registerState", None, [refparam(database, "dataBase"),
                                             refparam(state, "state")],
                     is_virtual=True)
        x.add_method("registerDerivatives", None, [refparam(database, "dataBase"),
                                                   refparam(derivatives, "derivatives")],
                     is_virtual=True)
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("enforceBoundaries", None, [refparam(state, "state"), refparam(derivatives, "derivatives")], is_virtual=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)
        
        # Attributes.
        const_ref_return_value(x, me, "%s::DdeviatoricStressDt" % me, symtensorfieldlist, [], "DdeviatoricStressDt")
        const_ref_return_value(x, me, "%s::bulkModulus" % me, scalarfieldlist, [], "bulkModulus")
        const_ref_return_value(x, me, "%s::shearModulus" % me, scalarfieldlist, [], "shearModulus")
        const_ref_return_value(x, me, "%s::yieldStrength" % me, scalarfieldlist, [], "yieldStrength")
        const_ref_return_value(x, me, "%s::plasticStrain0" % me, scalarfieldlist, [], "plasticStrain0")

        return
