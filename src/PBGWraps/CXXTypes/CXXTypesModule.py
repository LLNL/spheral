from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# Helper method to add std::vectors.
#-------------------------------------------------------------------------------
def generateStdVectorBindings(v, value, cppname, 
                              indexAsPointer = False,
                              wrapIterators = False):
    pointerValue = (value[-1:] == "*")

    v.add_constructor([])
    v.add_constructor([param("int", "size")])
    v.add_constructor([constrefparam(cppname, "rhs")])
    if not pointerValue:
        v.add_constructor([param("int", "size"), param(value, "value")])
    v.add_method("size", "int", [])
    v.add_method("size", "int", [], custom_name = "__len__")
    v.add_method("resize", None, [param("int", "size")])

    # Concatentation as addition.
    v.add_binary_numeric_operator("+")

    # Optionally expose the iterators.
    if wrapIterators:
        v.add_method("begin", "%s_iterator" % cppname, [])
        v.add_method("end", "%s_iterator" % cppname, [])

    if pointerValue:
        v.add_function_as_method("indexVectorOfPointers",
                                 retval(value, reference_existing_object=True),
                                 [param(cppname, "self"), param("int", "index")],
                                 template_parameters = [value],
                                 custom_name = "__getitem__")
    else:
        if indexAsPointer:
#             v.add_function_as_method("indexVectorAsReference",
#                                      retval(ref(value), reference_existing_object=True),
#                                      [param(cppname, "self"), param("int", "index")],
#                                      template_parameters = [value],
#                                      custom_name = "__getitem__")
            v.add_function_as_method("indexVectorAsPointer",
                                     retval(ptr(value), reference_existing_object=True),
                                     [param(cppname, "self"), param("int", "index")],
                                     template_parameters = [value],
                                     custom_name = "__getitem__")
        else:
            v.add_function_as_method("indexVector",
                                     value,
                                     [param(cppname, "self"), param("int", "index")],
                                     template_parameters = [value],
                                     custom_name = "__getitem__")

    if not pointerValue:
        v.add_function_as_method("assignToPosition", None, [param(cppname, "self"),
                                                            param("int", "index"),
                                                            param(value, "value")],
                                 template_parameters = [cppname],
                                 custom_name = "__setitem__")
        v.add_method("push_back", None, [param(value, "value")])
        v.add_method("push_back", None, [param(value, "value")], custom_name="append")
    else:
        v.add_function_as_method("appendToVectorOfPointers", None, [param(cppname, "self"),
                                                                    param(value, "value", transfer_ownership=False)],
                                 template_parameters = [cppname],
                                 custom_name = "append")

    return

#-------------------------------------------------------------------------------
# Helper method to add std::pairs.
#-------------------------------------------------------------------------------
def generateStdPairBindings(p, value1, value2, cppname,
                            extract_first = True,
                            extract_second = True):
    p.add_constructor([])
    val1ptr = value1[-1] == "*"
    val2ptr = value2[-1] == "*"

    # Constructors.
    p.add_constructor([])

    if val1ptr:
        firstparam = param(value1, "first", transfer_ownership=False)
    else:
        firstparam = param(value1, "first")

    if val2ptr:
        secondparam = param(value2, "second", transfer_ownership=False)
    else:
        secondparam = param(value2, "second")

    p.add_constructor([firstparam, secondparam])

    # First value.
    if extract_first:
        if val1ptr:
            p.add_function_as_method("extractFirstPairValue",
                                     retval(value1, reference_existing_object=True),
                                     [param(cppname, "self")],
                                     template_parameters = [cppname],
                                     custom_name = "first")
        else:
            p.add_instance_attribute("first", value1)

    # Second value.
    if extract_second:
        if val2ptr:
            p.add_function_as_method("extractSecondPairValue",
                                     retval(value1, reference_existing_object=True),
                                     [param(cppname, "self")],
                                     template_parameters = [cppname],
                                     custom_name = "second")
        else:
            p.add_instance_attribute("second", value2)

    return

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class CXXTypes:

    #---------------------------------------------------------------------------
    # Add all out stuff.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes
        mod.add_include('"CXXTypes/CXXTypes.hh"')

        # Namespace.
        std = mod.add_cpp_namespace("std")

        # These are the basic types we form vectors of.
        self.vectorElementTypes0 = ["char", "unsigned", "int", "float", "double", "string", "ULL"]
        self.vectorElementTypes1 = ["pair_unsigned_unsigned", "pair_ULL_ULL",
                                    "Vector1d", "Vector2d", "Vector3d",
                                    "Tensor1d", "Tensor2d", "Tensor3d",
                                    "SymTensor1d", "SymTensor2d", "SymTensor3d",
                                    "ThirdRankTensor1d", "ThirdRankTensor2d", "ThirdRankTensor3d",
                                    "Geom3Vector",
                                    "vector_of_char", "vector_of_unsigned", "vector_of_int", "vector_of_float", "vector_of_double", "vector_of_string", 
                                    "vector_of_Vector1d", "vector_of_Vector2d", "vector_of_Vector3d",
                                    "vector_of_Tensor1d", "vector_of_Tensor2d", "vector_of_Tensor3d",
                                    "vector_of_SymTensor1d", "vector_of_SymTensor2d", "vector_of_SymTensor3d",
                                    "vector_of_ThirdRankTensor1d", "vector_of_ThirdRankTensor2d", "vector_of_ThirdRankTensor3d",
                                    "vector_of_vector_of_char", "vector_of_vector_of_unsigned", "vector_of_vector_of_int", "vector_of_vector_of_float", "vector_of_vector_of_double", "vector_of_vector_of_string"]

        # The pair types.
        self.pairElementTypes = [("double", "double"), ("double", "string"), ("unsigned", "unsigned"), ("ULL", "ULL"), ("string", "string"),
                                 ("Vector1d", "Vector1d"), ("Vector2d", "Vector2d"), ("Vector3d", "Vector3d"), 
                                 ("Tensor1d", "Tensor1d"), ("Tensor2d", "Tensor2d"), ("Tensor3d", "Tensor3d")]

        # Vector types.
        for name in (self.vectorElementTypes0 + self.vectorElementTypes1):
            exec('self.vector_of_%s = addObject(mod, "vector_of_%s", allow_subclassing=True)' % (name, name))
            exec('self.vector_of_%s_iterator = addObject(mod, "vector_of_%s_iterator", allow_subclassing=True)' % (name, name))

        # Pair types.
        for (name1, name2) in self.pairElementTypes:
            exec('self.pair_%s_%s = addObject(mod, "pair_%s_%s", allow_subclassing=True)' % (name1, name2,
                                                                                             name1, name2))

        # A few value types need special mangling for the element specs.
        self.valueMap = {}
        for x in (self.vectorElementTypes0 + self.vectorElementTypes1):
            exec('self.valueMap["%s"] = "%s"' % (x, x))
        self.valueMap["unsigned"] = "unsigned int"
        self.valueMap["string"] = "std::string"
        self.valueMap["ULL"] = "uint64_t"

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        # Vector types.
        for name in self.vectorElementTypes0:
            exec('generateStdVectorBindings(self.vector_of_%s, "%s", "vector_of_%s")' % (name, self.valueMap[name], name))
        for name in self.vectorElementTypes1:
            exec('generateStdVectorBindings(self.vector_of_%s, "%s", "vector_of_%s", indexAsPointer=True)' % (name, name, name))

        # Pair types.
        for (name1, name2) in self.pairElementTypes:
            exec('generateStdPairBindings(self.pair_%s_%s, "%s", "%s", "pair_%s_%s")' % (name1, name2,
                                                                                         self.valueMap[name1], self.valueMap[name2],
                                                                                         name1, name2))
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["std"]

