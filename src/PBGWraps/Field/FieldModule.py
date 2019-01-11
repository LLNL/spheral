from pybindgen import *

from PBGutils import *
from CXXTypesModule import generateStdVectorBindings
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Field:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/FieldTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        self.FieldStorageType = space.add_enum("FieldStorageType", [("ReferenceFields", "Spheral::FieldStorageType::ReferenceFields"),
                                                                    ("CopyFields", "Spheral::FieldStorageType::CopyFields")])

        for ndim in self.dims:
            exec("""
# Expose types.
self.FieldBase%(dim)s = addObject(space, "FieldBase%(dim)s", allow_subclassing=True)
self.FieldListBase%(dim)s = addObject(space, "FieldListBase%(dim)s", allow_subclassing=True)

self.IntField%(dim)s =             addObject(space, "IntField%(dim)s", parent=self.FieldBase%(dim)s)
self.ULLField%(dim)s =             addObject(space, "ULLField%(dim)s", parent=self.FieldBase%(dim)s)
self.ScalarField%(dim)s =          addObject(space, "ScalarField%(dim)s", parent=self.FieldBase%(dim)s)
self.VectorField%(dim)s =          addObject(space, "VectorField%(dim)s", parent=self.FieldBase%(dim)s)
self.Vector3dField%(dim)s =        addObject(space, "Vector3dField%(dim)s", parent=self.FieldBase%(dim)s)
self.TensorField%(dim)s =          addObject(space, "TensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.SymTensorField%(dim)s =       addObject(space, "SymTensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.ThirdRankTensorField%(dim)s = addObject(space, "ThirdRankTensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.FourthRankTensorField%(dim)s = addObject(space, "FourthRankTensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.FifthRankTensorField%(dim)s = addObject(space, "FifthRankTensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.FacetedVolumeField%(dim)s =   addObject(space, "FacetedVolumeField%(dim)s", parent=self.FieldBase%(dim)s)
self.VectorDoubleField%(dim)s =    addObject(space, "VectorDoubleField%(dim)s", parent=self.FieldBase%(dim)s)
self.VectorVectorField%(dim)s =    addObject(space, "VectorVectorField%(dim)s", parent=self.FieldBase%(dim)s)
self.VectorTensorField%(dim)s =    addObject(space, "VectorTensorField%(dim)s", parent=self.FieldBase%(dim)s)
self.VectorSymTensorField%(dim)s = addObject(space, "VectorSymTensorField%(dim)s", parent=self.FieldBase%(dim)s)

self.IntFieldList%(dim)s =             addObject(space, "IntFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.ULLFieldList%(dim)s =             addObject(space, "ULLFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.ScalarFieldList%(dim)s =          addObject(space, "ScalarFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.VectorFieldList%(dim)s =          addObject(space, "VectorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.Vector3dFieldList%(dim)s =        addObject(space, "Vector3dFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.TensorFieldList%(dim)s =          addObject(space, "TensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.SymTensorFieldList%(dim)s =       addObject(space, "SymTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.ThirdRankTensorFieldList%(dim)s = addObject(space, "ThirdRankTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.FourthRankTensorFieldList%(dim)s = addObject(space, "FourthRankTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.FifthRankTensorFieldList%(dim)s = addObject(space, "FifthRankTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.FacetedVolumeFieldList%(dim)s =   addObject(space, "FacetedVolumeFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.VectorDoubleFieldList%(dim)s =    addObject(space, "VectorDoubleFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.VectorVectorFieldList%(dim)s =    addObject(space, "VectorVectorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.VectorTensorFieldList%(dim)s =    addObject(space, "VectorTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)
self.VectorSymTensorFieldList%(dim)s = addObject(space, "VectorSymTensorFieldList%(dim)s", parent=self.FieldListBase%(dim)s)

self.FieldListSet%(dim)s = addObject(space, "FieldListSet%(dim)s", allow_subclassing=True)
""" % {"dim" : "%id" % ndim})

            # std::vector<Field*>
            for element in ("Int", "Scalar", "Vector", "Tensor", "SymTensor"):
                exec("""
self.vector_of_%(element)sFieldPtr%(dim)s = addObject(mod, "vector_of_%(element)sFieldPtr%(dim)s", allow_subclassing=True)
self.vector_of_%(element)sFieldList%(dim)s = addObject(mod, "vector_of_%(element)sFieldList%(dim)s", allow_subclassing=True)
""" % {"element" : element, "dim" : "%id" % ndim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for ndim in self.dims:
            dim = "%id" % ndim
            polyvol = {1: "Box1d", 
                       2: "Polygon",
                       3: "Polyhedron"}[ndim]

            exec("""
self.addFieldBaseMethods(self.FieldBase%(dim)s, %(ndim)i)

self.addFieldMethods(self.IntField%(dim)s, "int", "IntField%(dim)s", %(ndim)i, applyNumberMethods=True)
self.addFieldMethods(self.ULLField%(dim)s, "unsigned long long", "ULLField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True)
self.addFieldMethods(self.ScalarField%(dim)s, "double", "ScalarField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True)
self.addFieldMethods(self.VectorField%(dim)s, "Vector%(dim)s", "VectorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.Vector3dField%(dim)s, "Geom3Vector", "Vector3dField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.TensorField%(dim)s, "Tensor%(dim)s", "TensorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.SymTensorField%(dim)s, "SymTensor%(dim)s", "SymTensorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.ThirdRankTensorField%(dim)s, "ThirdRankTensor%(dim)s", "ThirdRankTensorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.FourthRankTensorField%(dim)s, "FourthRankTensor%(dim)s", "FourthRankTensorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.FifthRankTensorField%(dim)s, "FifthRankTensor%(dim)s", "FifthRankTensorField%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.FacetedVolumeField%(dim)s, polyvol, "FacetedVolumeField%(dim)s", %(ndim)i, indexAsPointer=True)
self.addFieldMethods(self.VectorDoubleField%(dim)s, "vector_of_double", "VectorDoubleField%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.VectorVectorField%(dim)s, "vector_of_Vector%(dim)s", "VectorVectorField%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.VectorTensorField%(dim)s, "vector_of_Tensor%(dim)s", "VectorTensorField%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldMethods(self.VectorSymTensorField%(dim)s, "vector_of_SymTensor%(dim)s", "VectorSymTensorField%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)

self.addFieldListMethods(self.IntFieldList%(dim)s, "int", "IntFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True)
self.addFieldListMethods(self.ULLFieldList%(dim)s, "unsigned long long", "ULLFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True)
self.addFieldListMethods(self.ScalarFieldList%(dim)s, "double", "ScalarFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True)
self.addFieldListMethods(self.VectorFieldList%(dim)s, "Vector%(dim)s", "VectorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.Vector3dFieldList%(dim)s, "Geom3Vector", "Vector3dFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.TensorFieldList%(dim)s, "Tensor%(dim)s", "TensorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.SymTensorFieldList%(dim)s, "SymTensor%(dim)s", "SymTensorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.ThirdRankTensorFieldList%(dim)s, "ThirdRankTensor%(dim)s", "ThirdRankTensorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.FourthRankTensorFieldList%(dim)s, "FourthRankTensor%(dim)s", "FourthRankTensorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.FifthRankTensorFieldList%(dim)s, "FifthRankTensor%(dim)s", "FifthRankTensorFieldList%(dim)s", %(ndim)i, applyNumberMethods=True, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.FacetedVolumeFieldList%(dim)s, polyvol, "FacetedVolumeFieldList%(dim)s", %(ndim)i, indexAsPointer=True)
self.addFieldListMethods(self.VectorDoubleFieldList%(dim)s, "vector_of_double", "VectorDoubleFieldList%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.VectorVectorFieldList%(dim)s, "vector_of_Vector%(dim)s", "VectorVectorFieldList%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.VectorTensorFieldList%(dim)s, "vector_of_Tensor%(dim)s", "VectorTensorFieldList%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)
self.addFieldListMethods(self.VectorSymTensorFieldList%(dim)s, "vector_of_SymTensor%(dim)s", "VectorSymTensorFieldList%(dim)s", %(ndim)i, applyOrderingMethods=True, indexAsPointer=True)

self.addFieldListSetMethods(self.FieldListSet%(dim)s, %(ndim)i)
""" % {"dim" : dim, "ndim" : ndim})

            # std::vector<Field*> and std::vector<FieldList>
            for element in ["Int", "Scalar", "Vector", "Tensor", "SymTensor"]:
                exec("""
generateStdVectorBindings(self.vector_of_%(element)sFieldPtr%(dim)s, "Spheral::%(element)sField%(dim)s*", "vector_of_%(element)sFieldPtr%(dim)s", indexAsPointer=True)
generateStdVectorBindings(self.vector_of_%(element)sFieldList%(dim)s, "Spheral::%(element)sFieldList%(dim)s", "vector_of_%(element)sFieldList%(dim)s", indexAsPointer=True)
""" % {"element" : element, "dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add methods (Fields).
    #---------------------------------------------------------------------------
    def addFieldBaseMethods(self, x, ndim):

        # Object names.
        me = "Spheral::FieldBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
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
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim

        # # Constructors.
        # x.add_constructor([param("std::string", "name")])
        # x.add_constructor([param("std::string", "name"),
        #                    constrefparam(nodelist, "nodeList")])

        # Methods.
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
        x.add_method("size", "unsigned int", [], is_const=True, is_pure_virtual=True)
        x.add_method("Zero", None, [], is_pure_virtual=True)
        x.add_method("setNodeList", None, [constrefparam(nodelist, "nodeList")], is_pure_virtual=True)
        x.add_method("resizeField", None, [param("unsigned int", "size")], is_pure_virtual=True)
        x.add_method("resizeFieldInternal", None, [param("unsigned int", "size"),
                                                   param("unsigned int", "oldFirstGhostNode")], is_pure_virtual=True)
        x.add_method("resizeFieldGhost", None, [param("unsigned int", "size")], is_pure_virtual=True)
        x.add_method("deleteElement", None, [param("int", "nodeID")], is_pure_virtual=True)
        x.add_method("deleteElements", None, [constrefparam("vector_of_int", "nodeIDs")], is_pure_virtual=True)
        x.add_method("packValues", "vector_of_char", [constrefparam("vector_of_int", "nodeIDs")], is_const=True, is_pure_virtual=True)
        x.add_method("unpackValues", None, [param("int", "numElements"),
                                            param("int", "beginInsertionIndex"),
                                            constrefparam("vector_of_char", "buffer")], is_pure_virtual=True)

        # Attributes.
        x.add_instance_attribute("name", "std::string", getter="name", setter="name")

        return

    #---------------------------------------------------------------------------
    # Add methods (Fields).
    #---------------------------------------------------------------------------
    def addFieldMethods(self, x, val, me, ndim, applyNumberMethods=False, applyOrderingMethods=False, indexAsPointer=False):

        # Object names.
        me = "Spheral::%s" % me
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
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
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        if val == "unsigned long long":
            vector_of_value = "vector_of_ULL"
        else:
            vector_of_value = "vector_of_%s" % val

        # Constructors.
        x.add_constructor([param("std::string", "name")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(me, "field")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(nodelist, "nodeList")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(nodelist, "nodeList"),
                           param(val, "value")])
        x.add_constructor([constrefparam(me, "field")])

        # Methods.
        x.add_method("Zero", None, [], is_virtual=True)
        x.add_method("valid", "bool", [], is_const=True)
        x.add_method("internalValues", vector_of_value, [], is_const=True)
        x.add_method("ghostValues", vector_of_value, [], is_const=True)
        x.add_method("allValues", vector_of_value, [], is_const=True)

        # Comparison operators.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        if applyOrderingMethods:
            x.add_binary_comparison_operator("<")
            x.add_binary_comparison_operator(">")
            x.add_binary_comparison_operator("<=")
            x.add_binary_comparison_operator(">=")

        # These methods can only be used with Field values that can be used as numbers.
        if applyNumberMethods:
            x.add_method("applyMin", None, [param(val, "dataMin")])
            x.add_method("applyMax", None, [param(val, "dataMax")])
            x.add_method("sumElements", val, [], is_const=True)
            x.add_method("min", val, [], is_const=True)
            x.add_method("max", val, [], is_const=True)
            x.add_method("localSumElements", val, [], is_const=True)
            x.add_method("localMin", val, [], is_const=True)
            x.add_method("localMax", val, [], is_const=True)
            x.add_method("string", "std::string", [param("int", "precision", default_value="20")], is_const=True)
            x.add_method("string", None, [constrefparam("std::string", "string")])

            # We allow a few special methods for Symmetric tensors.
            if "SymTensor" in val:
                x.add_method("applyScalarMin", None, [param("double", "dataMin")])
                x.add_method("applyScalarMax", None, [param("double", "dataMax")])

        # Math operators.
        # These don't work yet, mostly 'cause the pybindgen implementation assumes that
        # the classes will be default constructable.
#         x.add_binary_numeric_operator("+")
#         x.add_binary_numeric_operator("-")
#         x.add_inplace_numeric_operator("+=")
#         x.add_inplace_numeric_operator("-=")
#         x.add_binary_numeric_operator("+", right_cppclass=val)
#         x.add_binary_numeric_operator("-", right_cppclass=val)
#         x.add_inplace_numeric_operator("+=", right_cppclass=val)
#         x.add_inplace_numeric_operator("-=", right_cppclass=val)

        # Sequence methods.
        x.add_method("size", "unsigned int", [], is_const=True, custom_name="__len__")
        if indexAsPointer:
            x.add_function_as_method("indexContainerAsPointer",
                                     retval(ptr(val), reference_existing_object=True),
                                     [param(me, "self"), param("int", "index")],
                                     template_parameters = [me],
                                     foreign_cpp_namespace = "Spheral",
                                     custom_name = "__getitem__")
        else:
            x.add_function_as_method("indexContainer", val,
                                     [param(me, "self"), param("int", "index")],
                                     template_parameters = [me],
                                     foreign_cpp_namespace = "Spheral",
                                     custom_name = "__getitem__")
        x.add_function_as_method("assignToContainerIndex", "int",
                                 [param(me, "self"), param("int", "index"), param(val, "value")],
                                 template_parameters = [me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "__setitem__")
        x.add_function_as_method("containsValue", "int",
                                 [param(me, "self"), param(val, "value")],
                                 template_parameters = [me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "__contains__")

        # The virtual methods from FieldBase.
        x.add_method("size", "unsigned int", [], is_const=True, is_virtual=True)
        x.add_method("Zero", None, [], is_virtual=True)
        x.add_method("setNodeList", None, [constrefparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("resizeField", None, [param("unsigned int", "size")], is_virtual=True)
        x.add_method("resizeFieldInternal", None, [param("unsigned int", "size"),
                                                   param("unsigned int", "oldFirstGhostNode")], is_virtual=True)
        x.add_method("resizeFieldGhost", None, [param("unsigned int", "size")], is_virtual=True)
        x.add_method("deleteElement", None, [param("int", "nodeID")], is_virtual=True)
        x.add_method("deleteElements", None, [constrefparam("vector_of_int", "nodeIDs")], is_virtual=True)
        x.add_method("packValues", "vector_of_char", [constrefparam("vector_of_int", "nodeIDs")], is_const=True, is_virtual=True)
        x.add_method("unpackValues", None, [param("int", "numElements"),
                                            param("int", "beginInsertionIndex"),
                                            constrefparam("vector_of_char", "buffer")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("numElements", "int", is_const=True, getter="numElements")
        x.add_instance_attribute("numInternalElements", "int", is_const=True, getter="numInternalElements")

        return

    #---------------------------------------------------------------------------
    # Add methods (FieldLists).
    #---------------------------------------------------------------------------
    def addFieldListMethods(self, x, val, me, ndim, applyNumberMethods=False, applyOrderingMethods=False, indexAsPointer=False):

        # Object names.
        me = "Spheral::%s" % me
        field = me.replace("List", "")
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldBase%id" % ndim
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
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        vector_of_nodelist = "vector_of_NodeList%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("FieldStorageType", "aStorageType")])
        x.add_constructor([constrefparam(me, "rhs")])

        # Methods.
        x.add_method("copyFields", None, [])
        x.add_method("haveField", "bool", [constrefparam(field, "field")], is_const=True)
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")], is_const=True)
        x.add_method("assignFields", None, [constrefparam(me, "fieldList")])
        x.add_method("appendField", None, [constrefparam(field, "field")])
        x.add_method("deleteField", None, [constrefparam(field, "field")])
        x.add_method("appendNewField", None, [param("std::string", "name"), constrefparam(nodelist, "nodeList"), param(val, "value")])
        x.add_method("setMasterNodeLists", None, [constrefparam(vector, "r"), constrefparam(symtensor, "H"),
                                                  refparam("vector_of_vector_of_int", "masterLists"),
                                                  refparam("vector_of_vector_of_int", "coarseNeighbors")], is_const=True)
        x.add_method("setMasterNodeLists", None, [constrefparam(vector, "r"),
                                                  refparam("vector_of_vector_of_int", "masterLists"),
                                                  refparam("vector_of_vector_of_int", "coarseNeighbors")], is_const=True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "r"), constrefparam(symtensor, "H"),
                                                  constrefparam("vector_of_vector_of_int", "coarseNeighbors"),
                                                  refparam("vector_of_vector_of_int", "refineNeighbors")], is_const=True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "r"),
                                                  constrefparam("vector_of_vector_of_int", "coarseNeighbors"),
                                                  refparam("vector_of_vector_of_int", "refineNeighbors")], is_const=True)
        x.add_method("Zero", None, [])
        x.add_method("nodeListPtrs", vector_of_nodelist, [], is_const=True)

        # Comparison operators.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        if applyOrderingMethods:
            x.add_binary_comparison_operator("<")
            x.add_binary_comparison_operator(">")
            x.add_binary_comparison_operator("<=")
            x.add_binary_comparison_operator(">=")

        # These methods can only be used with Field values that can be used as numbers.
        if applyNumberMethods:
            x.add_method("applyMin", None, [param(val, "dataMin")])
            x.add_method("applyMax", None, [param(val, "dataMax")])
            x.add_method("sumElements", val, [], is_const=True)
            x.add_method("min", val, [], is_const=True)
            x.add_method("max", val, [], is_const=True)
            x.add_method("localSumElements", val, [], is_const=True)
            x.add_method("localMin", val, [], is_const=True)
            x.add_method("localMax", val, [], is_const=True)

            # We allow a few special methods for Symmetric tensors.
            if "SymTensor" in val:
                x.add_method("applyScalarMin", None, [param("double", "dataMin")])
                x.add_method("applyScalarMax", None, [param("double", "dataMax")])

        # Sequence methods.
        x.add_method("size", "unsigned int", [], is_const=True, custom_name="__len__")
        x.add_function_as_method("indexContainer",
                                 retval(ptr(field), reference_existing_object=True),
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "__getitem__")
        x.add_function_as_method("assignToContainerIndexPtr", "int",
                                 [param(me, "self"), param("int", "index"), param(ptr(field), "value", transfer_ownership=False)],
                                 template_parameters = [me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "__setitem__")
        x.add_function_as_method("containsValue", "int",
                                 [param(me, "self"), param(ptr(field), "value", transfer_ownership=False)],
                                 template_parameters = [me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "__contains__")

        x.add_function_as_method("fieldForNodeList",
                                 retval(ptr(field), reference_existing_object=True),
                                 [param(me, "self"), 
                                  constrefparam(nodelist, "nodeList")],
                                 template_parameters = [dim, me],
                                 custom_name = "fieldForNodeList")

        # Extract an individual value.
        if indexAsPointer:
            x.add_function_as_method("indexFieldListForValuePointer",
                                     retval(ptr(val), reference_existing_object=True),
                                     [param(me, "self"), param("int", "fieldIndex"), param("int", "nodeIndex")],
                                     template_parameters = [me],
                                     custom_name = "__call__")
        else:
            x.add_function_as_method("indexFieldListForValue",
                                     val,
                                     [param(me, "self"), param("int", "fieldIndex"), param("int", "nodeIndex")],
                                     template_parameters = [me],
                                     custom_name = "__call__")

        # Attributes.
        x.add_instance_attribute("storageType", "FieldStorageType", is_const=True, getter="storageType")
        x.add_instance_attribute("numFields", "int", is_const=True, getter="numFields")
        x.add_instance_attribute("numNodes", "int", is_const=True, getter="numNodes")
        x.add_instance_attribute("numInternalNodes", "int", is_const=True, getter="numInternalNodes")
        x.add_instance_attribute("numGhostNodes", "int", is_const=True, getter="numGhostNodes")

        return

    #---------------------------------------------------------------------------
    # Add methods (FieldListSet).
    #---------------------------------------------------------------------------
    def addFieldListSetMethods(self, x, ndim):

        # Object names.
        me = "FieldListSet%id" % ndim
        vector_of_scalarfieldlist = "vector_of_ScalarFieldList%id" % ndim
        vector_of_vectorfieldlist = "vector_of_VectorFieldList%id" % ndim
        vector_of_tensorfieldlist = "vector_of_TensorFieldList%id" % ndim
        vector_of_symtensorfieldlist = "vector_of_SymTensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])

        # Attributes.
        x.add_instance_attribute("ScalarFieldLists", retval(ptr(vector_of_scalarfieldlist), reference_existing_object=True), getter="ScalarFieldListPtrs", is_const=True)
        x.add_instance_attribute("VectorFieldLists", retval(ptr(vector_of_vectorfieldlist), reference_existing_object=True), getter="VectorFieldListPtrs", is_const=True)
        x.add_instance_attribute("TensorFieldLists", retval(ptr(vector_of_tensorfieldlist), reference_existing_object=True), getter="TensorFieldListPtrs", is_const=True)
        x.add_instance_attribute("SymTensorFieldLists", retval(ptr(vector_of_symtensorfieldlist), reference_existing_object=True), getter="SymTensorFieldListPtrs", is_const=True)

        return

