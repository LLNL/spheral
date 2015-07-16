from pybindgen import *

from PBGutils import *
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class DataBase:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/DataBaseTypes.hh"' % srcdir)
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("DataBaseSpace")

        # Expose types.
        self.StateBase1d = addObject(Spheral, "StateBase1d", allow_subclassing=True)
        self.StateBase2d = addObject(Spheral, "StateBase2d", allow_subclassing=True)
        self.StateBase3d = addObject(Spheral, "StateBase3d", allow_subclassing=True)

        self.State1d = addObject(Spheral, "State1d", allow_subclassing=True, parent=self.StateBase1d)
        self.State2d = addObject(Spheral, "State2d", allow_subclassing=True, parent=self.StateBase2d)
        self.State3d = addObject(Spheral, "State3d", allow_subclassing=True, parent=self.StateBase3d)
        
        self.StateDerivatives1d = addObject(Spheral, "StateDerivatives1d", allow_subclassing=True, parent=self.StateBase1d)
        self.StateDerivatives2d = addObject(Spheral, "StateDerivatives2d", allow_subclassing=True, parent=self.StateBase2d)
        self.StateDerivatives3d = addObject(Spheral, "StateDerivatives3d", allow_subclassing=True, parent=self.StateBase3d)
        
        self.DataBase1d = addObject(space, "DataBase1d", allow_subclassing=True)
        self.DataBase2d = addObject(space, "DataBase2d", allow_subclassing=True)
        self.DataBase3d = addObject(space, "DataBase3d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.addStateBaseMethods(self.StateBase1d, 1)
        self.addStateBaseMethods(self.StateBase2d, 2)
        self.addStateBaseMethods(self.StateBase3d, 3)

        self.addStateMethods(self.State1d, 1)
        self.addStateMethods(self.State2d, 2)
        self.addStateMethods(self.State3d, 3)

        self.addStateDerivativesMethods(self.StateDerivatives1d, 1)
        self.addStateDerivativesMethods(self.StateDerivatives2d, 2)
        self.addStateDerivativesMethods(self.StateDerivatives3d, 3)

        self.addDataBaseMethods(self.DataBase1d, 1)
        self.addDataBaseMethods(self.DataBase2d, 2)
        self.addDataBaseMethods(self.DataBase3d, 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["DataBaseSpace"]

    #---------------------------------------------------------------------------
    # Add StateBase methods.
    #---------------------------------------------------------------------------
    def addStateBaseMethods(self, x, ndim):

        # Object names.
        me = "StateBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        thirdranktensor = "ThirdRankTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        fieldlistbase = "Spheral::FieldSpace::FieldListBase%id" % ndim
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
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        mesh = ("Spheral::MeshSpace::LineMesh",
                "Spheral::MeshSpace::PolygonalMesh",
                "Spheral::MeshSpace::PolyhedralMesh")[ndim - 1]

        # Constructors.
        x.add_constructor([])

        # Comparison operators.
        x.add_binary_comparison_operator("==")

        # Methods.
        x.add_method("registered", "bool", [constrefparam("std::string", "key")], is_const=True)
        x.add_method("registered", "bool", [constrefparam(fieldbase, "field")], is_const=True)
        x.add_method("fieldNameRegistered", "bool", [constrefparam("std::string", "fieldName")], is_const=True)
        x.add_method("enroll", None, [refparam(fieldbase, "field")], is_virtual=True)
        x.add_method("enroll", None, [refparam(fieldlistbase, "fieldList")], is_virtual=True)

        # Extract Fields.
        for result, value, customname in [(intfield, "int", "intField"),
                                          (scalarfield, "double", "scalarField"),
                                          (vectorfield, vector, "vectorField"),
                                          (tensorfield, tensor, "tensorField"),
                                          (symtensorfield, symtensor, "symTensorField")]:
            exec("""x.add_function_as_method("fieldFromStateBase",
                                 retval(ptr("%(result)s"), reference_existing_object=True),
                                 [param(me, "self"), constrefparam("std::string", "key")],
                                 template_parameters = [dim, "%(value)s"],
                                 custom_name = "%(customname)s")""" %
                 {"result" : result, "value" : value, "customname" : customname})


        # Extract FieldLists.
        for value, default, result, customname in [("int", "0", intfieldlist, "intFields"),
                                                   ("double", "0.0", scalarfieldlist, "scalarFields"),
                                                   (vector, "%s::zero" % vector, vectorfieldlist, "vectorFields"),
                                                   (tensor, "%s::zero" % tensor, tensorfieldlist, "tensorFields"),
                                                   (symtensor, "%s::zero" % symtensor, symtensorfieldlist, "symTensorFields"),
                                                   (thirdranktensor, "%s::zero" % thirdranktensor, thirdranktensorfieldlist, "thirdRankTensorFields"),
                                                   ("vector_of_double", "vector_of_double()", vectordoublefieldlist, "vectorDoubleFields"),
                                                   ("vector_of_%s" % vector, "vector_of_%s()" % vector, vectorvectorfieldlist, "vectorVectorFields"),
                                                   ("vector_of_%s" % symtensor, "vector_of_%s()" % symtensor, vectorsymtensorfieldlist, "vectorSymTensorFields")]:
            exec("""x.add_method("fields", "%(result)s", 
                     [param("std::string", "name"), param("%(value)s", "dummy", default_value="%(default)s")],
                     template_parameters = ["%(value)s"], 
                     custom_name = "%(customname)s",
                     is_const = True)""" %
                 {"value" : value, "default" : default, "result" : result, "customname" : customname})

        # Extract allFields, set of fields for a given data type.
        for value, default, result, customname in [("int", "0", "vector_of_IntFieldPtr%id" % ndim, "allIntFields"),
                                                   ("double", "0.0", "vector_of_ScalarFieldPtr%id" % ndim, "allScalarFields"),
                                                   (vector, "%s::zero" % vector, "vector_of_VectorFieldPtr%id" % ndim, "allVectorFields"),
                                                   (tensor, "%s::zero" % tensor, "vector_of_TensorFieldPtr%id" % ndim, "allTensorFields"),
                                                   (symtensor, "%s::zero" % symtensor, "vector_of_SymTensorFieldPtr%id" % ndim, "allSymTensorFields")]:
            exec("""x.add_method("allFields", "%(result)s",
                     [param("%(value)s", "dummy", default_value="%(default)s")],
                     template_parameters = ["%(value)s"],
                     custom_name = "%(customname)s",
                     is_const = True)""" %
                 {"value" : value, "default" : default, "result" : result, "customname" : customname})

        x.add_method("keys", "vector_of_string", [], is_const=True)
        x.add_method("fieldKeys", "vector_of_string", [], is_const=True)

        x.add_method("meshRegistered", "bool", [], is_const=True)
        x.add_method("mesh", mesh, [])

        x.add_method("assign", None, [constrefparam(me, "rhs")])
        x.add_method("copyState", None, [], is_virtual=True)

        x.add_method("key", "std::string", [param(fieldbase, "field")], is_static=True)
        x.add_method("buildFieldKey", "std::string", 
                     [param("std::string", "fieldName"), param("std::string", "nodeListName")],
                     is_static = True)
        x.add_method("splitFieldKey", None,
                     [param("std::string", "key"),
                      refparam("std::string", "fieldKey"),
                      refparam("std::string", "nodeListKey")],
                     is_static = True)

        return

    #---------------------------------------------------------------------------
    # Add State methods.
    #---------------------------------------------------------------------------
    def addStateMethods(self, x, ndim):

        # Object names.
        me = "State%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        vectorpkgs = "vector_of_Physics%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(database, "dataBase"),
                           refparam(vectorpkgs, "packages")])
        x.add_constructor([constrefparam(me, "state")])

        # Comparison operators.
        x.add_binary_comparison_operator("==")

        # Methods.
        x.add_method("update", None, [refparam(derivatives, "derivatives"),
                                      param("double", "multiplier"),
                                      param("double", "t"),
                                      param("double", "dt")])
        x.add_method("policyKeys", "vector_of_string", [], is_const=True)
        
        # Attributes.
        x.add_instance_attribute("timeAdvanceOnly", "bool", getter="timeAdvanceOnly", setter="timeAdvanceOnly")

        return

    #---------------------------------------------------------------------------
    # Add StateDerivatives methods.
    #---------------------------------------------------------------------------
    def addStateDerivativesMethods(self, x, ndim):

        # Object names.
        me = "StateDerivatives%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        vectorpkgs = "vector_of_Physics%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(database, "dataBase"),
                           refparam(vectorpkgs, "packages")])
        x.add_constructor([constrefparam(me, "rhs")])

        # Comparison operators.
        x.add_binary_comparison_operator("==")

        # Methods.
        x.add_method("initializeNodePairInformation", None, [])
        x.add_method("calculatedNodePairsSymmetric", "bool", [], is_const=True)
        x.add_method("Zero", None, [])

        return

    #---------------------------------------------------------------------------
    # Add methods.
    #---------------------------------------------------------------------------
    def addDataBaseMethods(self, x, ndim):

        # Object names.
        me = "Spheral::DataBaseSpace::DataBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        polyvol = {1: "Box1d", 
                   2: "Polygon",
                   3: "Polyhedron"}[ndim]
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        polyvolfield = "Spheral::FieldSpace::FacetedVolumeField%id" % ndim
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
        polyvolfieldlist = "Spheral::FieldSpace::FacetedVolumeFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        fluidnodelist = "Spheral::NodeSpace::FluidNodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        vectorpkgs = "vector_of_Physics%id" % ndim
        vectornodelists = "vector_of_NodeList%id" % ndim
        vectorfluidnodelists = "vector_of_FluidNodeList%id" % ndim
        vector_of_Vector = "vector_of_Vector%id" % ndim

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("updateConnectivityMap", None, [param("bool", "computeGhostConnectivity")], is_const=True)
        x.add_method("patchConnectivityMap", None, [constrefparam(intfieldlist, "flags"),
                                                    constrefparam(intfieldlist, "old2new")], is_const=True)
        const_ref_return_value(x, me, "%s::connectivityMap" % me, connectivitymap, [], "connectivityMap")
        x.add_function_as_method("connectivityMapFromDataBase", 
                                 retval(ptr(connectivitymap), reference_existing_object=True),
                                 [param(me, "self"), param("bool", "computeGhostConnectivity")],
                                 template_parameters = [dim],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "connectivityMap")
        x.add_method("appendNodeList", None, [refparam(fluidnodelist, "nodeList")])
        x.add_method("appendNodeList", None, [refparam(nodelist, "nodeList")])
        x.add_method("deleteNodeList", None, [refparam(fluidnodelist, "nodeList")])
        x.add_method("deleteNodeList", None, [refparam(nodelist, "nodeList")])
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")])

        x.add_method("nodeListPtrs", vectornodelists, [], is_const=True, custom_name="nodeLists")
        x.add_method("fluidNodeListPtrs", vectorfluidnodelists, [], is_const=True, custom_name="fluidNodeLists")

##         x.add_function_as_method("const_reference_as_pointer",
##                                  retval(ptr(vectornodelists), reference_existing_object=True),
##                                  [param(me, "self")],
##                                  template_parameters = [me, vectornodelists, "&%s::nodeListPtrs" % me],
##                                  foreign_cpp_namespace = "Spheral",
##                                  custom_name = "nodeLists")
##         x.add_function_as_method("const_reference_as_pointer",
##                                  retval(ptr(vectorfluidnodelists), reference_existing_object=True),
##                                  [param(me, "self")],
##                                  template_parameters = [me, vectorfluidnodelists, "&%s::fluidNodeListPtrs" % me],
##                                  foreign_cpp_namespace = "Spheral",
##                                  custom_name = "fluidNodeLists")

        x.add_method("setMasterNodeLists", None, [constrefparam(vector, "position"),
                                                  constrefparam(symtensor, "H")],
                     is_const = True)
        x.add_method("setMasterFluidNodeLists", None, [constrefparam(vector, "position"),
                                                       constrefparam(symtensor, "H")],
                     is_const = True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "position"),
                                                  constrefparam(symtensor, "H")],
                     is_const = True)
        x.add_method("setRefineFluidNodeLists", None, [constrefparam(vector, "position"),
                                                       constrefparam(symtensor, "H")],
                     is_const = True)

        for result, value, default, customname in [(intfieldlist, "int", "0", "Int"),
                                                   (scalarfieldlist, "double", "0.0", "Scalar"),
                                                   (vectorfieldlist, vector, "%s::zero" % vector, "Vector"),
                                                   (tensorfieldlist, tensor, "%s::zero" % tensor, "Tensor"),
                                                   (symtensorfieldlist, symtensor, "%s::zero" % symtensor, "SymTensor"),
                                                   (polyvolfieldlist, polyvol, "%s()" % polyvol, "FacetedVolume")]:
            exec("""
x.add_method("newGlobalFieldList", "%(result)s", [param("%(value)s", "value", default_value="%(default)s"),
                                                  param("std::string", "name", default_value='"unnamed field list"')],
             template_parameters = ["%(value)s"],
             is_const = True,
             custom_name = "newGlobal%(customname)sFieldList")
x.add_method("newFluidFieldList", "%(result)s", [param("%(value)s", "value", default_value="%(default)s"),
                                                 param("std::string", "name", default_value='"unnamed field list"')],
             template_parameters = ["%(value)s"],
             is_const = True,
             custom_name = "newFluid%(customname)sFieldList")
""" % {"result" : result, "value" : value, "default" : default, "customname" : customname})

        x.add_method("boundingBox", None, [refparam(vector, "xmin"),
                                           refparam(vector, "xmax"),
                                           param("bool", "ghost", default_value="true")], is_const=True)
        x.add_method("boundingBox", None, [refparam(vector, "xmin"),
                                           refparam(vector, "xmax"),
                                           constrefparam(intfieldlist, "mask"),
                                           param("bool", "ghost", default_value="true")], is_const=True)
        x.add_method("localSamplingBoundingVolume", None, [refparam(vector, "centroid"),
                                                           refparam("double", "radiusNodes"),
                                                           refparam("double", "radiusSample"),
                                                           refparam(vector, "xminNodes"),
                                                           refparam(vector, "xmaxNodes"),
                                                           refparam(vector, "xminSample"),
                                                           refparam(vector, "xmaxSample")], is_const=True)
        x.add_method("globalSamplingBoundingVolume", None, [refparam(vector, "centroid"),
                                                            refparam("double", "radiusNodes"),
                                                            refparam("double", "radiusSample"),
                                                            refparam(vector, "xminNodes"),
                                                            refparam(vector, "xmaxNodes"),
                                                            refparam(vector, "xminSample"),
                                                            refparam(vector, "xmaxSample")], is_const=True)
        x.add_method("localSamplingBoundingBoxes", None, [refparam(vector_of_Vector, "xminima"), refparam(vector_of_Vector, "xmaxima")], is_const=True)
        x.add_method("globalSamplingBoundingBoxes", None, [refparam(vector_of_Vector, "xminima"), refparam(vector_of_Vector, "xmaxima")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True)

        x.add_method("globalHinverse", None, [refparam(symtensorfieldlist, "result")], is_const=True)
        x.add_method("fluidHinverse", None, [refparam(symtensorfieldlist, "result")], is_const=True)
        x.add_method("fluidPressure", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidTemperature", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidSoundSpeed", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidVolume", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidLinearMomentum", None, [refparam(vectorfieldlist, "result")], is_const=True)
        x.add_method("fluidTotalEnergy", None, [refparam(scalarfieldlist, "result")], is_const=True)

        # Attributes.
        x.add_instance_attribute("numNodeLists", "int", getter="numNodeLists", is_const=True)
        x.add_instance_attribute("numFluidNodeLists", "int", getter="numFluidNodeLists", is_const=True)
        x.add_instance_attribute("numInternalNodes", "int", getter="numInternalNodes", is_const=True)
        x.add_instance_attribute("numGhostNodes", "int", getter="numGhostNodes", is_const=True)
        x.add_instance_attribute("numNodes", "int", getter="numNodes", is_const=True)
        x.add_instance_attribute("globalNumInternalNodes", "int", getter="globalNumInternalNodes", is_const=True)
        x.add_instance_attribute("globalNumGhostNodes", "int", getter="globalNumGhostNodes", is_const=True)
        x.add_instance_attribute("globalNumNodes", "int", getter="globalNumNodes", is_const=True)

        x.add_instance_attribute("globalMass", scalarfieldlist, getter="globalMass", is_const=True)
        x.add_instance_attribute("globalPosition", vectorfieldlist, getter="globalPosition", is_const=True)
        x.add_instance_attribute("globalVelocity", vectorfieldlist, getter="globalVelocity", is_const=True)
        x.add_instance_attribute("globalHfield", symtensorfieldlist, getter="globalHfield", is_const=True)
        x.add_instance_attribute("globalWork", scalarfieldlist, getter="globalWork", is_const=True)

        x.add_instance_attribute("fluidMass", scalarfieldlist, getter="fluidMass", is_const=True)
        x.add_instance_attribute("fluidPosition", vectorfieldlist, getter="fluidPosition", is_const=True)
        x.add_instance_attribute("fluidVelocity", vectorfieldlist, getter="fluidVelocity", is_const=True)
        x.add_instance_attribute("fluidMassDensity", scalarfieldlist, getter="fluidMassDensity", is_const=True)
        x.add_instance_attribute("fluidSpecificThermalEnergy", scalarfieldlist, getter="fluidSpecificThermalEnergy", is_const=True)
        x.add_instance_attribute("fluidHfield", symtensorfieldlist, getter="fluidHfield", is_const=True)
        x.add_instance_attribute("fluidWork", scalarfieldlist, getter="fluidWork", is_const=True)

        x.add_instance_attribute("globalNodeExtent", vectorfieldlist, getter="globalNodeExtent", is_const=True)
        x.add_instance_attribute("fluidNodeExtent", vectorfieldlist, getter="fluidNodeExtent", is_const=True)

        x.add_instance_attribute("numNeighbors", intfieldlist, getter="numNeighbors", is_const=True)

        x.add_instance_attribute("maxKernelExtent", "double", getter="maxKernelExtent", is_const=True)

        x.add_instance_attribute("nDim", "int")
        x.add_static_attribute("nDim", "int", True)

        return
