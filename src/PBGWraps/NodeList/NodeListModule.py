from pybindgen import *

from PBGutils import *
from ref_return_value import *

import sys
sys.path.append("../CXXTypes")
from CXXTypesModule import generateStdPairBindings, generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class NodeList:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"NodeList/NodeListTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral.add_cpp_namespace("NodeSpace")

        # Expose types.
        self.NodeType = self.space.add_enum("NodeType", ["InternalNode", "GhostNode"])

        for dim, ndim in (("1d", 1), 
                          ("2d", 2),
                          ("3d", 3)):
            exec("""
self.NodeListRegistrar%(dim)s = addObject(Spheral, "NodeListRegistrar%(dim)s", is_singleton=True);
self.NodeList%(dim)s = addObject(self.space, "NodeList%(dim)s", allow_subclassing=True)
self.FluidNodeList%(dim)s = addObject(self.space, "FluidNodeList%(dim)s", allow_subclassing=True, parent=self.NodeList%(dim)s)
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

""" % {"dim" : dim, "ndim" : ndim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        Spheral = mod.add_cpp_namespace("Spheral")

        for dim, ndim, mesh in (("1d", 1, "LineMesh"),
                                ("2d", 2, "PolygonalMesh"),
                                ("3d", 3, "PolyhedralMesh")):
            exec("""
self.addNodeListRegistrarMethods(self.NodeListRegistrar%(dim)s, %(ndim)s)
self.addNodeListMethods(self.NodeList%(dim)s, %(ndim)s)
self.addFluidNodeListMethods(self.FluidNodeList%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseMethods(self.SmoothingScaleBase%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.FixedSmoothingScale%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.SPHSmoothingScale%(dim)s, %(ndim)s)
self.addSmoothingScaleBaseDescendentMethods(self.ASPHSmoothingScale%(dim)s, %(ndim)s)

generateStdVectorBindings(self.vector_of_NodeList%(dim)s, "Spheral::NodeSpace::NodeList%(dim)s*", "vector_of_NodeList%(dim)s")
generateStdVectorBindings(self.vector_of_FluidNodeList%(dim)s, "Spheral::NodeSpace::FluidNodeList%(dim)s*", "vector_of_FluidNodeList%(dim)s")

generateStdPairBindings(self.pair_NodeList%(dim)s_string,
                        "const Spheral::NodeSpace::NodeList%(dim)s*",
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
                         constrefparam(mesh, "mesh"),
                         constrefparam("Vector%(dim)s", "xmin"),
                         constrefparam("Vector%(dim)s", "xmax"),
                         param("unsigned int", "numInternal"),
                         param("double", "nPerh"),
                         param("double", "voidThreshold"),
                         refparam("Spheral::NodeSpace::NodeList%(dim)s", "voidNodes")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>"],
                        custom_name = "generateVoidNodes",
                        docstring = "Generate a new set of void nodes based on surfaces in a set of NodeLists.")

self.space.add_function("nthNodalMoment", "Spheral::FieldSpace::ScalarFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::KernelSpace::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "0U"],
                        custom_name = "zerothNodalMoment",
                        docstring = "Compute the zeroth moment of the local node distribution in eta space.")
self.space.add_function("nthNodalMoment", "Spheral::FieldSpace::VectorFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::KernelSpace::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "1U"],
                        custom_name = "firstNodalMoment",
                        docstring = "Compute the first moment of the local node distribution in eta space.")
self.space.add_function("nthNodalMoment", "Spheral::FieldSpace::SymTensorFieldList%(dim)s",
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::KernelSpace::TableKernel%(dim)s", "xmin"),
                         param("bool", "renormalize")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>", "2U"],
                        custom_name = "secondNodalMoment",
                        docstring = "Compute the second moment of the local node distribution in eta space.")
self.space.add_function("zerothAndFirstNodalMoments", None,
                        [constrefparam("vector_of_NodeList%(dim)s", "nodeLists"),
                         constrefparam("Spheral::KernelSpace::TableKernel%(dim)s", "xmin"),
                         param("bool", "useGradientAsKernel"),
                         refparam("Spheral::FieldSpace::ScalarFieldList%(dim)s", "zerothMoment"),
                         refparam("Spheral::FieldSpace::VectorFieldList%(dim)s", "firstMoment")],
                        template_parameters = ["Spheral::Dim<%(ndim)s>"],
                        custom_name = "zerothAndFirstNodalMoments",
                        docstring = "Compute the zeroth and first moments of the local node distribution in eta space.")

""" % {"dim" : dim, "ndim" : ndim})

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
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        neighbor = "Spheral::NeighborSpace::Neighbor%id" % ndim

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

        me = "Spheral::NodeSpace::NodeList%id" % ndim

        # External objects.
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        neighbor = "Spheral::NeighborSpace::Neighbor%id" % ndim

        # Constructors.
        x.add_constructor([param("std::string", "name"),
                           param("int", "numInternal", default_value="0"),
                           param("int", "numGhost", default_value="0"),
                           param("double", "hmin", default_value="1.0e-20"),
                           param("double", "hmax", default_value="1.0e20"),
                           param("double", "hminratio", default_value="0.1"),
                           param("double", "nPerh", default_value="2.01"),
                           param("int", "maxNumNeighbors", default_value="500")])

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
        x.add_instance_attribute("numNodes", "int", getter="numNodes", is_const=True)
        x.add_instance_attribute("numInternalNodes", "int", getter="numInternalNodes", setter="numInternalNodes")
        x.add_instance_attribute("numGhostNodes", "int", getter="numGhostNodes", setter="numGhostNodes")
        x.add_instance_attribute("numFields", "int", getter="numFields", is_const=True)
        x.add_instance_attribute("firstGhostNode", "int", getter="firstGhostNode", is_const=True)
        x.add_instance_attribute("nodesPerSmoothingScale", "double", getter="nodesPerSmoothingScale", setter="nodesPerSmoothingScale")
        x.add_instance_attribute("maxNumNeighbors", "int", getter="maxNumNeighbors", setter="maxNumNeighbors")
        x.add_instance_attribute("hmin", "double", getter="hmin", setter="hmin")
        x.add_instance_attribute("hmax", "double", getter="hmax", setter="hmax")
        x.add_instance_attribute("hminratio", "double", getter="hminratio", setter="hminratio")

        return

    #---------------------------------------------------------------------------
    # Add methods to FluidNodeList.
    #---------------------------------------------------------------------------
    def addFluidNodeListMethods(self, x, ndim):

        me = "Spheral::NodeSpace::FluidNodeList%id" % ndim

        # External objects.
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        fluidderivativeproducer = "Spheral::NodeSpace::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::Material::EquationOfState%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

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
    # Add methods to SmoothingScaleBase.
    #---------------------------------------------------------------------------
    def addSmoothingScaleBaseMethods(self, x, ndim):

        me = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        # External objects.
        vector = "Spheral::Vector%id" % ndim
        tensor = "Spheral::Tensor%id" % ndim
        symtensor = "Spheral::SymTensor%id" % ndim
        fluidnodelist = "Spheral::NodeSpace::FluidNodeList%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        fluidderivativeproducer = "Spheral::NodeSpace::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::Material::EquationOfState%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([])
        #x.add_constructor([constrefparam(me, "rhs")])

        # Methods.
        x.add_method("newSmoothingScaleAndDerivative", None, [constrefparam(symtensorfield, "H"),
                                                              constrefparam(tensorfield, "DvDx"),
                                                              constrefparam(scalarfield, "zerothMoment"),
                                                              constrefparam(symtensorfield, "secondMoment"),
                                                              constrefparam(connectivitymap, "connectivityMap"),
                                                              constrefparam(tablekernel, "W"),
                                                              param("double", "hmin"),
                                                              param("double", "hmax"),
                                                              param("double", "hminratio"),
                                                              param("double", "nPerh"),
                                                              param("int", "maxNumNeighbors"),
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
        fluidnodelist = "Spheral::NodeSpace::FluidNodeList%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        fluidderivativeproducer = "Spheral::NodeSpace::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::Material::EquationOfState%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

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
        fluidnodelist = "Spheral::NodeSpace::FluidNodeList%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        state = "Spheral::State%id" % ndim
        statederivatives = "Spheral::StateDerivatives%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"
        fluidderivativeproducer = "Spheral::NodeSpace::FluidDerivativeProducer%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
        equationofstate = "Spheral::Material::EquationOfState%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("smoothingScaleDerivative", symtensor, [constrefparam(symtensor, "H"),
                                                             constrefparam(tensor, "DvDx"),
                                                             param("double", "hmin"),
                                                             param("double", "hmax"),
                                                             param("double", "hminratio"),
                                                             param("double", "nPerh"),
                                                             param("int", "maxNumNeighbors")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("newSmoothingScale", symtensor, [constrefparam(symtensor, "H"),
                                                      param("double", "zerothMoment"),
                                                      constrefparam(symtensor, "secondMoment"),
                                                      param("int", "numNeighbors"),
                                                      constrefparam(tablekernel, "W"),
                                                      param("double", "hmin"),
                                                      param("double", "hmax"),
                                                      param("double", "hminratio"),
                                                      param("double", "nPerh"),
                                                      param("int", "maxNumNeighbors")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
        x.add_method("idealSmoothingScale", symtensor, [constrefparam(symtensor, "H"),
                                                        param("double", "zerothMoment"),
                                                        constrefparam(symtensor, "secondMoment"),
                                                        param("int", "numNeighbors"),
                                                        constrefparam(tablekernel, "W"),
                                                        param("double", "hmin"),
                                                        param("double", "hmax"),
                                                        param("double", "hminratio"),
                                                        param("double", "nPerh"),
                                                        param("int", "maxNumNeighbors")],
                     is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

        return

