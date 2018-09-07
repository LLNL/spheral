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
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/DataBaseTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        for dim in self.dims:
            exec('''
self.StateBase%(dim)id = addObject(space, "StateBase%(dim)id", allow_subclassing=True)
self.State%(dim)id = addObject(space, "State%(dim)id", allow_subclassing=True, parent=self.StateBase%(dim)id)
self.StateDerivatives%(dim)id = addObject(space, "StateDerivatives%(dim)id", allow_subclassing=True, parent=self.StateBase%(dim)id)
self.DataBase%(dim)id = addObject(space, "DataBase%(dim)id", allow_subclassing=True)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for dim in self.dims:
            exec('''
self.addStateBaseMethods(self.StateBase%(dim)id, %(dim)i)
self.addStateMethods(self.State%(dim)id, %(dim)i)
self.addStateDerivativesMethods(self.StateDerivatives%(dim)id, %(dim)i)
self.addDataBaseMethods(self.DataBase%(dim)id, %(dim)i)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

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
        fieldbase = "Spheral::FieldBase%id" % ndim
        fieldlistbase = "Spheral::FieldListBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        boundary = "Spheral::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        mesh = ("Spheral::LineMesh",
                "Spheral::PolygonalMesh",
                "Spheral::PolyhedralMesh")[ndim - 1]

        # Constructors.
        x.add_constructor([])

        # Comparison operators.
        x.add_binary_comparison_operator("==")

        # Methods.
        x.add_method("registered", "bool", [param("std::string", "key")], is_const=True)
        x.add_method("registered", "bool", [param(fieldbase, "field")], is_const=True)
        x.add_method("fieldNameRegistered", "bool", [param("std::string", "fieldName")], is_const=True)
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
        database = "Spheral::DataBase%id" % ndim
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
        database = "Spheral::DataBase%id" % ndim
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
    # Add DataBase methods.
    #---------------------------------------------------------------------------
    def addDataBaseMethods(self, x, ndim):

        # Object names.
        me = "Spheral::DataBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        thirdranktensor = "ThirdRankTensor%id" % ndim
        fourthranktensor = "FourthRankTensor%id" % ndim
        fifthranktensor = "FifthRankTensor%id" % ndim
        polyvol = {1: "Box1d", 
                   2: "Polygon",
                   3: "Polyhedron"}[ndim]
        fieldbase = "Spheral::FieldBase%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        vector3dfield = "Spheral::Vector3dField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        fourthranktensorfield = "Spheral::FourthRankTensorField%id" % ndim
        fifthranktensorfield = "Spheral::FifthRankTensorField%id" % ndim
        polyvolfield = "Spheral::FacetedVolumeField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        fourthranktensorfieldlist = "Spheral::FourthRankTensorFieldList%id" % ndim
        fifthranktensorfieldlist = "Spheral::FifthRankTensorFieldList%id" % ndim
        polyvolfieldlist = "Spheral::FacetedVolumeFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        fluidnodelist = "Spheral::FluidNodeList%id" % ndim
        solidnodelist = "Spheral::SolidNodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        vectorpkgs = "vector_of_Physics%id" % ndim
        vectornodelists = "vector_of_NodeList%id" % ndim
        vectorfluidnodelists = "vector_of_FluidNodeList%id" % ndim
        vectorsolidnodelists = "vector_of_SolidNodeList%id" % ndim
        vector_of_Vector = "vector_of_Vector%id" % ndim

        # Constructors.
        x.add_constructor([])

        # Methods.
        x.add_method("reinitializeNeighbors", None, [], is_const=True)
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
        x.add_method("appendNodeList", None, [refparam(solidnodelist, "nodeList")])
        x.add_method("appendNodeList", None, [refparam(fluidnodelist, "nodeList")])
        x.add_method("appendNodeList", None, [refparam(nodelist, "nodeList")])
        x.add_method("deleteNodeList", None, [refparam(solidnodelist, "nodeList")])
        x.add_method("deleteNodeList", None, [refparam(fluidnodelist, "nodeList")])
        x.add_method("deleteNodeList", None, [refparam(nodelist, "nodeList")])
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")])

        x.add_method("nodeListPtrs", vectornodelists, [], is_const=True, custom_name="nodeLists")
        x.add_method("fluidNodeListPtrs", vectorfluidnodelists, [], is_const=True, custom_name="fluidNodeLists")
        x.add_method("solidNodeListPtrs", vectorsolidnodelists, [], is_const=True, custom_name="solidNodeLists")

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
                                                  constrefparam(symtensor, "H"),
                                                  refparam("vector_of_vector_of_int", "masterLists"),
                                                  refparam("vector_of_vector_of_int", "coarseNeighbors")],
                     is_const = True)
        x.add_method("setMasterFluidNodeLists", None, [constrefparam(vector, "position"),
                                                       constrefparam(symtensor, "H"),
                                                       refparam("vector_of_vector_of_int", "masterLists"),
                                                       refparam("vector_of_vector_of_int", "coarseNeighbors")],
                     is_const = True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "position"),
                                                  constrefparam(symtensor, "H"),
                                                  constrefparam("vector_of_vector_of_int", "coarseNeighbors"),
                                                  refparam("vector_of_vector_of_int", "refineNeighbors")],
                     is_const = True)
        x.add_method("setRefineFluidNodeLists", None, [constrefparam(vector, "position"),
                                                       constrefparam(symtensor, "H"),
                                                       constrefparam("vector_of_vector_of_int", "coarseNeighbors"),
                                                       refparam("vector_of_vector_of_int", "refineNeighbors")],
                     is_const = True)

        for result, value, default, customname in [(intfieldlist, "int", "0", "Int"),
                                                   (scalarfieldlist, "double", "0.0", "Scalar"),
                                                   (vectorfieldlist, vector, "%s::zero" % vector, "Vector"),
                                                   (tensorfieldlist, tensor, "%s::zero" % tensor, "Tensor"),
                                                   (symtensorfieldlist, symtensor, "%s::zero" % symtensor, "SymTensor"),
                                                   (thirdranktensorfieldlist, thirdranktensor, "%s::zero" % thirdranktensor, "ThirdRankTensor"),
                                                   (fourthranktensorfieldlist, fourthranktensor, "%s::zero" % fourthranktensor, "FourthRankTensor"),
                                                   (fifthranktensorfieldlist, fifthranktensor, "%s::zero" % fifthranktensor, "FifthRankTensor"),
                                                   (polyvolfieldlist, polyvol, "%s()" % polyvol, "FacetedVolume"),
                                                   (vectordoublefieldlist, "vector_of_double", "vector_of_double()", "vector_of_double"),
                                                   (vectorvectorfieldlist, vector_of_Vector, vector_of_Vector + "()", "vector_of_Vector")]:
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
x.add_method("newSolidFieldList", "%(result)s", [param("%(value)s", "value", default_value="%(default)s"),
                                                 param("std::string", "name", default_value='"unnamed field list"')],
             template_parameters = ["%(value)s"],
             is_const = True,
             custom_name = "newSolid%(customname)sFieldList")
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
        x.add_method("fluidGamma", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidEntropy", None, [refparam(scalarfieldlist, "result")], is_const=True)
        x.add_method("fluidLinearMomentum", None, [refparam(vectorfieldlist, "result")], is_const=True)
        x.add_method("fluidTotalEnergy", None, [refparam(scalarfieldlist, "result")], is_const=True)

        # Attributes.
        x.add_instance_attribute("numNodeLists", "int", getter="numNodeLists", is_const=True)
        x.add_instance_attribute("numFluidNodeLists", "int", getter="numFluidNodeLists", is_const=True)
        x.add_instance_attribute("numSolidNodeLists", "int", getter="numSolidNodeLists", is_const=True)
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

        x.add_instance_attribute("solidDeviatoricStress", symtensorfieldlist, getter="solidDeviatoricStress", is_const=True)
        x.add_instance_attribute("solidPlasticStrain", scalarfieldlist, getter="solidPlasticStrain", is_const=True)
        x.add_instance_attribute("solidPlasticStrainRate", scalarfieldlist, getter="solidPlasticStrainRate", is_const=True)
        x.add_instance_attribute("solidDamage", symtensorfieldlist, getter="solidDamage", is_const=True)
        x.add_instance_attribute("solidEffectiveDamage", symtensorfieldlist, getter="solidEffectiveDamage", is_const=True)
        x.add_instance_attribute("solidDamageGradient", vectorfieldlist, getter="solidDamageGradient", is_const=True)
        x.add_instance_attribute("solidFragmentIDs", intfieldlist, getter="solidFragmentIDs", is_const=True)

        x.add_instance_attribute("globalNodeExtent", vectorfieldlist, getter="globalNodeExtent", is_const=True)
        x.add_instance_attribute("fluidNodeExtent", vectorfieldlist, getter="fluidNodeExtent", is_const=True)

        x.add_instance_attribute("numNeighbors", intfieldlist, getter="numNeighbors", is_const=True)

        x.add_instance_attribute("maxKernelExtent", "double", getter="maxKernelExtent", is_const=True)

        x.add_instance_attribute("nDim", "int", is_const=True)
        x.add_instance_attribute("isRZ", "bool", is_const=True)

        x.add_static_attribute("nDim", "int", True)
        x.add_static_attribute("isRZ", "bool", True)

        return
