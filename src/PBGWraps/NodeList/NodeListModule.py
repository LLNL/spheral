from pybindgen import *

from PBGutils import *
from ref_return_value import *

from CXXTypesModule import generateStdPairBindings, generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class NodeList:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/NodeListTypes.hh"' % srcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        self.NodeType = self.space.add_enum("NodeType", [("InternalNode", "Spheral::NodeType::InternalNode"), 
                                                         ("GhostNode",    "Spheral::NodeType::GhostNode")])

        for ndim in self.dims:
            exec("""
self.NodeListRegistrar%(dim)s = addObject(self.space, "NodeListRegistrar%(dim)s", is_singleton=True);
self.NodeList%(dim)s = addObject(self.space, "NodeList%(dim)s", allow_subclassing=True)
self.FluidNodeList%(dim)s = addObject(self.space, "FluidNodeList%(dim)s", allow_subclassing=True, parent=self.NodeList%(dim)s)
self.SolidNodeList%(dim)s = addObject(self.space, "SolidNodeList%(dim)s", allow_subclassing=True, parent=self.FluidNodeList%(dim)s)
self.SmoothingScaleBase%(dim)s = addObject(self.space, "SmoothingScaleBase%(dim)s", allow_subclassing=True)

self.FixedSmoothingScale%(dim)s = addObject(self.space, "FixedSmoothingScale%(dim)s", allow_subclassing=True, parent=self.SmoothingScaleBase%(dim)s)
self.SPHSmoothingScale%(dim)s = addObject(self.space, "SPHSmoothingScale%(dim)s", allow_subclassing=True, parent=self.SmoothingScaleBase%(dim)s)
self.ASPHSmoothingScale%(dim)s = addObject(self.space, "ASPHSmoothingScale%(dim)s", allow_subclassing=True, parent=self.SmoothingScaleBase%(dim)s)

self.vector_of_NodeList%(dim)s = addObject(mod, "vector_of_NodeList%(dim)s", allow_subclassing=True)
self.vector_of_FluidNodeList%(dim)s = addObject(mod, "vector_of_FluidNodeList%(dim)s", allow_subclassing=True)

self.pair_NodeList%(dim)s_string = addObject(mod, "pair_NodeList%(dim)s_string", allow_subclassing=True)
        
self.vector_of_pair_NodeList%(dim)s_string = addObject(mod, "vector_of_pair_NodeList%(dim)s_string", allow_subclassing=True)

self.vector_of_NodeList%(dim)s_iterator = addObject(mod, "vector_of_NodeList%(dim)s_iterator", allow_subclassing=True)
self.vector_of_FluidNodeList%(dim)s_iterator = addObject(mod, "vector_of_FluidNodeList%(dim)s_iterator", allow_subclassing=True)
self.vector_of_SolidNodeList%(dim)s = addObject(mod, "vector_of_SolidNodeList%(dim)s", allow_subclassing=True)
self.vector_of_SolidNodeList%(dim)s_iterator = addObject(mod, "vector_of_SolidNodeList%(dim)s_iterator", allow_subclassing=True)

""" % {"ndim" : ndim,
       "dim"  : "%id" % ndim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        Spheral = mod.add_cpp_namespace("Spheral")

        for ndim in self.dims:
            exec("""
self.addNodeListRegistrarMethods(self.NodeListRegistrar%(dim)s, %(ndim)s)
self.addNodeListMethods(self.NodeList%(dim)s, %(ndim)s)
self.addFluidNodeListMethods(self.FluidNodeList%(dim)s, %(ndim)s)
self.addSolidNodeListBindings(self.SolidNodeList%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseMethods(self.SmoothingScaleBase%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.FixedSmoothingScale%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.SPHSmoothingScale%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.ASPHSmoothingScale%(dim)s, %(ndim)s)

generateStdVectorBindings(self.vector_of_NodeList%(dim)s, "Spheral::NodeList%(dim)s*", "vector_of_NodeList%(dim)s")
generateStdVectorBindings(self.vector_of_FluidNodeList%(dim)s, "Spheral::FluidNodeList%(dim)s*", "vector_of_FluidNodeList%(dim)s")
generateStdVectorBindings(self.vector_of_SolidNodeList%(dim)s, "Spheral::SolidNodeList%(dim)s*", "vector_of_SolidNodeList%(dim)s")

generateStdPairBindings(self.pair_NodeList%(dim)s_string,
                        "const Spheral::NodeList%(dim)s*",
                        "std::string",
                        "pair_NodeList%(dim)s_string",
                        extract_first = False)

generateStdVectorBindings(self.vector_of_pair_NodeList%(dim)s_string, "pair_NodeList%(dim)s_string", "vector_of_pair_NodeList%(dim)s_string")

Spheral.add_function("nodeListRegistrarInstance",
                     retval(ptr("Spheral::NodeListRegistrar%(dim)s"), caller_owns_return=True),
                     [],
                     template_parameters = ["Spheral::Dim<%(ndim)s>"],
                     custom_name = "nodeListRegistrarInstance%(dim)s")

self.space.add_function("generateVoidNodes", None,
                        [constrefparam("vector_of_Vector%(dim)s", "generators"),
                         constrefparam("vector_of_SymTensor%(dim)s", "Hs"),
                         constrefparam("%(mesh)s", "mesh"),
                         constrefparam("Vector%(dim)s", "xmin"),
                         constrefparam("Vector%(dim)s", "xmax"),
                         param("unsigned int", "numInternal"),
                         param("double", "nPerh"),
                         param("double", "voidThreshold"),
                         refparam("Spheral::NodeList%(dim)s", "voidNodes")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>"],
                        custom_name = "generateVoidNodes",
                        docstring = "Generate a new set of void nodes based on surfaces in a set of NodeLists.")

self.space.add_function("nthNodalMoment", "Spheral::ScalarFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "0U"],
                        custom_name = "zerothNodalMoment",
                        docstring = "Compute the zeroth moment of the local node distribution in eta space.")
self.space.add_function("nthNodalMoment", "Spheral::VectorFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "1U"],
                        custom_name = "firstNodalMoment",
                        docstring = "Compute the first moment of the local node distribution in eta space.")
self.space.add_function("nthNodalMoment", "Spheral::SymTensorFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "2U"],
                        custom_name = "secondNodalMoment",
                        docstring = "Compute the second moment of the local node distribution in eta space.")
self.space.add_function("zerothAndFirstNodalMoments", None,
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::TableKernel%(dim)s", "xmin"),
                         param("bool", "useGradientAsKernel"),
                         refparam("Spheral::ScalarFieldList%(dim)s", "zerothMoment"),
                         refparam("Spheral::VectorFieldList%(dim)s", "firstMoment")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>"],
                        custom_name = "zerothAndFirstNodalMoments",
                        docstring = "Compute the zeroth and first moments of the local node distribution in eta space.")

""" % {"ndim" : ndim,
       "dim"  : "%id" % ndim,
       "mesh" : ("LineMesh", "PolygonalMesh", "PolyhedralMesh")[ndim - 1]})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["NodeSpace"]

    #---------------------------------------------------------------------------
    # Add methods to NodeListRegistrar.
    #---------------------------------------------------------------------------
    def addNodeListRegistrarMethods(self, x, ndim):

        # External objects.
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        fileio = "Spheral::FileIO"
        neighbor = "Spheral::Neighbor%id" % ndim

        # Methods.
        x.add_method("valid", "bool", [], is_const=True)

        # Attributes.
        x.add_instance_attribute("numNodeLists", "int", getter="numNodeLists", is_const=True)
        x.add_instance_attribute("numFluidNodeLists", "int", getter="numFluidNodeLists", is_const=True)
        x.add_instance_attribute("domainDecompositionIndependent", "bool", getter="domainDecompositionIndependent", setter="domainDecompositionIndependent")
        x.add_instance_attribute("registeredNames", "vector_of_string", getter="registeredNames", is_const=True)

        return

    #---------------------------------------------------------------------------
    # Add methods to NodeList.
    #---------------------------------------------------------------------------
    def addNodeListMethods(self, x, ndim):

        me = "Spheral::NodeList%id" % ndim

        # External objects.
        fieldbase = "Spheral::FieldBase%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        fileio = "Spheral::FileIO"
        neighbor = "Spheral::Neighbor%id" % ndim

        # Constructors.
        x.add_constructor([param("std::string", "name"),
                           param("unsigned int", "numInternal", default_value="0"),
                           param("unsigned int", "numGhost", default_value="0"),
                           param("double", "hmin", default_value="1.0e-20"),
                           param("double", "hmax", default_value="1.0e20"),
                           param("double", "hminratio", default_value="0.1"),
                           param("double", "nPerh", default_value="2.01"),
                           param("unsigned int", "maxNumNeighbors", default_value="500")])

        # Methods.
        x.add_method("nodeType", "NodeType", [param("int", "ID")], is_const=True)
        const_ref_return_value(x, me, "%s::mass" % me, scalarfield, [], "mass")
        const_ref_return_value(x, me, "%s::positions" % me, vectorfield, [], "positions")
        const_ref_return_value(x, me, "%s::velocity" % me, vectorfield, [], "velocity")
        const_ref_return_value(x, me, "%s::Hfield" % me, symtensorfield, [], "Hfield")
        ref_return_value(x, me, "%s::work" % me, scalarfield, [], "work")
        
        x.add_method("mass", None, [constrefparam(scalarfield, "newValue")])
        x.add_method("positions", None, [constrefparam(vectorfield, "newValue")])
        x.add_method("velocity", None, [constrefparam(vectorfield, "newValue")])
        x.add_method("Hfield", None, [constrefparam(symtensorfield, "newValue")])
        x.add_method("work", None, [constrefparam(scalarfield, "newValue")])

        x.add_method("Hinverse", None, [refparam(symtensorfield, "result")], is_const=True)

        ref_return_value(x, me, "%s::neighbor" % me, neighbor, [], "neighbor")
        x.add_method("registerNeighbor", None, [refparam(neighbor, "neighbor")])
        x.add_method("unregisterNeighbor", None, [])

        x.add_method("haveField", "bool", [constrefparam(fieldbase, "field")], is_const=True)

        # Virtual methods.
        x.add_method("deleteNodes", None, [constrefparam("vector_of_int", "nodeIDs")], is_virtual=True)
        x.add_method("reorderNodes", None, [constrefparam("vector_of_int", "newOrdering")], is_virtual=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("name", "std::string", getter="name", is_const=True)
        x.add_instance_attribute("numNodes", "unsigned int", getter="numNodes", is_const=True)
        x.add_instance_attribute("numInternalNodes", "unsigned int", getter="numInternalNodes", setter="numInternalNodes")
        x.add_instance_attribute("numGhostNodes", "unsigned int", getter="numGhostNodes", setter="numGhostNodes")
        x.add_instance_attribute("numFields", "unsigned int", getter="numFields", is_const=True)
        x.add_instance_attribute("firstGhostNode", "unsigned int", getter="firstGhostNode", is_const=True)
        x.add_instance_attribute("nodesPerSmoothingScale", "double", getter="nodesPerSmoothingScale", setter="nodesPerSmoothingScale")
        x.add_instance_attribute("maxNumNeighbors", "unsigned int", getter="maxNumNeighbors", setter="maxNumNeighbors")
        x.add_instance_attribute("hmin", "double", getter="hmin", setter="hmin")
        x.add_instance_attribute("hmax", "double", getter="hmax", setter="hmax")
        x.add_instance_attribute("hminratio", "double", getter="hminratio", setter="hminratio")

        # Comparison operators.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        return

    #---------------------------------------------------------------------------
    # Add methods to FluidNodeList.
    #---------------------------------------------------------------------------
    def addFluidNodeListMethods(self, x, ndim):

        me = "Spheral::FluidNodeList%id" % ndim

        # External objects.
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIO"
        fluidderivativeproducer = "Spheral::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::EquationOfState%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([param("std::string", "name"),
                           refparam(equationofstate, "eos"),
                           param("int", "numInternal", default_value="0"),
                           param("int", "numGhost", default_value="0"),
                           param("double", "hmin", default_value="1.0e-20"),
                           param("double", "hmax", default_value="1.0e20"),
                           param("double", "hminratio", default_value="0.1"),
                           param("double", "nPerh", default_value="2.01"),
                           param("int", "maxNumNeighbors", default_value="500"),
                           param("double", "rhoMin", default_value="1.0e-10"),
                           param("double", "rhoMax", default_value="1.0e100")])

        # Methods.
        const_ref_return_value(x, me, "%s::massDensity" % me, scalarfield, [], "massDensity")
        const_ref_return_value(x, me, "%s::specificThermalEnergy" % me, scalarfield, [], "specificThermalEnergy")

        x.add_method("massDensity", None, [constrefparam(scalarfield, "newValue")])
        x.add_method("specificThermalEnergy", None, [constrefparam(scalarfield, "newValue")])

        x.add_method("pressure", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("temperature", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("soundSpeed", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("volume", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("linearMomentum", None, [refparam(vectorfield, "result")], is_const=True, is_virtual=True)
        x.add_method("totalEnergy", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)

        const_ref_return_value(x, me, "%s::equationOfState" % me, equationofstate, [], "equationOfState")
        x.add_method("equationOfState", None, [constrefparam(equationofstate, "equationOfState")])

        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")],
                     is_const = True,
                     is_virtual = True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")],
                     is_virtual = True)

        # Attributes.
        x.add_instance_attribute("rhoMin", "double", getter="rhoMin", setter="rhoMin")
        x.add_instance_attribute("rhoMax", "double", getter="rhoMax", setter="rhoMax")

    #---------------------------------------------------------------------------
    # SolidNodeList
    #---------------------------------------------------------------------------
    def addSolidNodeListBindings(self, x, ndim):

        me = "Spheral::SolidNodeList%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::EquationOfState%id" % ndim
        strengthmodel = "Spheral::StrengthModel%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([param("std::string", "name"),
                           refparam(equationofstate, "eos"),
                           refparam(strengthmodel, "strength"),
                           param("int", "numInternal", default_value="0"),
                           param("int", "numGhost", default_value="0"),
                           param("double", "hmin", default_value="0.0"),
                           param("double", "hmax", default_value="1.0e100"),
                           param("double", "hminratio", default_value="0.1"),
                           param("double", "nPerh", default_value="2.01"),
                           param("int", "maxNumNeighbors", default_value="500"),
                           param("double", "rhoMin", default_value="1.0e-10"),
                           param("double", "rhoMax", default_value="1.0e100")])

        # Methods.
        x.add_method("soundSpeed", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("bulkModulus", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("shearModulus", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)
        x.add_method("yieldStrength", None, [refparam(scalarfield, "result")], is_const=True, is_virtual=True)

        const_ref_return_value(x, me, "%s::deviatoricStress" % me, symtensorfield, [], "deviatoricStress")
        const_ref_return_value(x, me, "%s::plasticStrain" % me, scalarfield, [], "plasticStrain")
        const_ref_return_value(x, me, "%s::plasticStrainRate" % me, scalarfield, [], "plasticStrainRate")
        const_ref_return_value(x, me, "%s::damage" % me, symtensorfield, [], "damage")
        const_ref_return_value(x, me, "%s::effectiveDamage" % me, symtensorfield, [], "effectiveDamage")
        const_ref_return_value(x, me, "%s::damageGradient" % me, vectorfield, [], "damageGradient")
        const_ref_return_value(x, me, "%s::fragmentIDs" % me, intfield, [], "fragmentIDs")
        const_ref_return_value(x, me, "%s::strengthModel" % me, strengthmodel, [], "strengthModel")

        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [constrefparam(fileio, "FileIO"), refparam("std::string", "pathName")], is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # Add methods to SmoothingScaleBase.
    #---------------------------------------------------------------------------
    def addSmoothingScaleBaseMethods(self, x, ndim):

        me = "Spheral::SmoothingScaleBase%id" % ndim

        # External objects.
        vector = "Spheral::Vector%id" % ndim
        tensor = "Spheral::Tensor%id" % ndim
        symtensor = "Spheral::SymTensor%id" % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIO"
        fluidderivativeproducer = "Spheral::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::EquationOfState%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([])
        #x.add_constructor([constrefparam(me, "rhs")])

        # Methods.
        x.add_method("newSmoothingScaleAndDerivative", None, [constrefparam(symtensorfield, "H"),
                                                              constrefparam(vectorfield, "position"),
                                                              constrefparam(tensorfield, "DvDx"),
                                                              constrefparam(scalarfield, "zerothMoment"),
                                                              constrefparam(symtensorfield, "secondMoment"),
                                                              constrefparam(connectivitymap, "connectivityMap"),
                                                              constrefparam(tablekernel, "W"),
                                                              param("double", "hmin"),
                                                              param("double", "hmax"),
                                                              param("double", "hminratio"),
                                                              param("double", "nPerh"),
                                                              refparam(symtensorfield, "DHDt"),
                                                              refparam(symtensorfield, "Hideal")], is_const=True)

        self.addSmoothingScaleBaseVirtualMethods(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # Add methods to SmoothingScaleBase descendents.
    #---------------------------------------------------------------------------
    def addSmoothingScaleBaseDescendentMethods(self, x, ndim):

        # External objects.
        vector = "Spheral::Vector%id" % ndim
        tensor = "Spheral::Tensor%id" % ndim
        symtensor = "Spheral::SymTensor%id" % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIO"
        fluidderivativeproducer = "Spheral::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::EquationOfState%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([])

        # Methods.
        self.addSmoothingScaleBaseVirtualMethods(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # Add virtual methods to SmoothingScaleBase types.
    #---------------------------------------------------------------------------
    def addSmoothingScaleBaseVirtualMethods(self, x, ndim, pureVirtual):

        # External objects.
        vector = "Spheral::Vector%id" % ndim
        tensor = "Spheral::Tensor%id" % ndim
        symtensor = "Spheral::SymTensor%id" % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldlist%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldlist%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIO"
        fluidderivativeproducer = "Spheral::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::EquationOfState%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"
        mesh = "Spheral::" + {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]
        zone = "%s::Zone" % mesh

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("smoothingScaleDerivative", symtensor, [constrefparam(symtensor, "H"),
                                                             constrefparam(vector, "pos"),
                                                             constrefparam(tensor, "DvDx"),
                                                             param("double", "hmin"),
                                                             param("double", "hmax"),
                                                             param("double", "hminratio"),
                                                             param("double", "nPerh")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("newSmoothingScale", symtensor, [constrefparam(symtensor, "H"),
                                                      constrefparam(vector, "pos"),
                                                      param("double", "zerothMoment"),
                                                      constrefparam(symtensor, "secondMoment"),
                                                      constrefparam(tablekernel, "W"),
                                                      param("double", "hmin"),
                                                      param("double", "hmax"),
                                                      param("double", "hminratio"),
                                                      param("double", "nPerh"),
                                                      constrefparam(connectivitymap, "connectivityMap"),
                                                      param("unsigned int", "nodeListi"),
                                                      param("unsigned int", "i")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("idealSmoothingScale", symtensor, [constrefparam(symtensor, "H"),
                                                        constrefparam(vector, "pos"),
                                                        param("double", "zerothMoment"),
                                                        constrefparam(symtensor, "secondMoment"),
                                                        constrefparam(tablekernel, "W"),
                                                        param("double", "hmin"),
                                                        param("double", "hmax"),
                                                        param("double", "hminratio"),
                                                        param("double", "nPerh"),
                                                        constrefparam(connectivitymap, "connectivityMap"),
                                                        param("unsigned int", "nodeListi"),
                                                        param("unsigned int", "i")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("idealSmoothingScale", symtensor, [constrefparam(symtensor, "H"),
                                                        constrefparam(mesh, "mesh"),
                                                        constrefparam(zone, "zone"),
                                                        param("double", "hmin"),
                                                        param("double", "hmax"),
                                                        param("double", "hminratio"),
                                                        param("double", "nPerh")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

        return

