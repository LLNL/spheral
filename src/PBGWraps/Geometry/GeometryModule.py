from pybindgen import *

from PBGutils import *
from CXXTypesModule import generateStdVectorBindings

double = Parameter.new("double", "double")

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Geometry:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/GeometryTypes.hh"' % srcdir)
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral

        # Expose types.
        for dim in ("1d", "2d", "3d"):
            exec("""
self.Vector%(dim)s = addObject(Spheral, "Vector%(dim)s")
self.Tensor%(dim)s = addObject(Spheral, "Tensor%(dim)s")
self.SymTensor%(dim)s = addObject(Spheral, "SymTensor%(dim)s")
self.ThirdRankTensor%(dim)s = addObject(Spheral, "ThirdRankTensor%(dim)s")
self.FourthRankTensor%(dim)s = addObject(Spheral, "FourthRankTensor%(dim)s")
self.FifthRankTensor%(dim)s = addObject(Spheral, "FifthRankTensor%(dim)s")
self.EigenStruct%(dim)s = addObject(Spheral, "EigenStruct%(dim)s")
self.Plane%(dim)s = addObject(Spheral, "Plane%(dim)s")
self.vector_of_FacetedVolume%(dim)s = addObject(mod, "vector_of_FacetedVolume%(dim)s", allow_subclassing=True)
self.vector_of_vector_of_FacetedVolume%(dim)s = addObject(mod, "vector_of_vector_of_FacetedVolume%(dim)s", allow_subclassing=True)
self.vector_of_Plane%(dim)s = addObject(mod, "vector_of_Plane%(dim)s", allow_subclassing=True)
""" % {"dim" : dim})

        self.Geom3Vector = addObject(Spheral, "Geom3Vector")
        
        self.Box1d = addObject(Spheral, "Box1d")

        self.Facet2d = addObject(Spheral, "Facet2d")
        self.Polygon = addObject(Spheral, "Polygon")

        self.Facet3d = addObject(Spheral, "Facet3d")
        self.Polyhedron = addObject(Spheral, "Polyhedron")

        self.vector_of_Facet2d = addObject(mod, "vector_of_Facet2d", allow_subclassing=True)
        self.vector_of_Facet3d = addObject(mod, "vector_of_Facet3d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        # Add methods to Vectors.
        self.addVectorMethods(self.Vector1d, 1)
        self.addVectorMethods(self.Vector2d, 2)
        self.addVectorMethods(self.Vector3d, 3)
        self.addGeom3VectorMethods(self.Geom3Vector)

        # Add methods to Tensors.
        self.Tensor1d.add_constructor([param("double", "xx", default_value = "0.0")])
        self.Tensor2d.add_constructor([param("double", "xx", default_value = "0.0"),
                                       param("double", "xy", default_value = "0.0"),
                                       param("double", "yx", default_value = "0.0"),
                                       param("double", "yy", default_value = "0.0")])
        self.Tensor3d.add_constructor([param("double", "xx", default_value = "0.0"),
                                       param("double", "xy", default_value = "0.0"),
                                       param("double", "xz", default_value = "0.0"),
                                       param("double", "yx", default_value = "0.0"),
                                       param("double", "yy", default_value = "0.0"),
                                       param("double", "yz", default_value = "0.0"),
                                       param("double", "zx", default_value = "0.0"),
                                       param("double", "zy", default_value = "0.0"),
                                       param("double", "zz", default_value = "0.0")])
        self.addTensorMethods(self.Tensor1d, "Tensor1d", 1)
        self.addTensorMethods(self.Tensor2d, "Tensor2d", 2)
        self.addTensorMethods(self.Tensor3d, "Tensor3d", 3)

        # Add methods to SymTensors.
        self.SymTensor1d.add_constructor([param("double", "xx", default_value = "0.0")])
        self.SymTensor2d.add_constructor([param("double", "xx", default_value = "0.0"),
                                          param("double", "xy", default_value = "0.0"),
                                          param("double", "yx", default_value = "0.0"),
                                          param("double", "yy", default_value = "0.0")])
        self.SymTensor3d.add_constructor([param("double", "xx", default_value = "0.0"),
                                          param("double", "xy", default_value = "0.0"),
                                          param("double", "xz", default_value = "0.0"),
                                          param("double", "yx", default_value = "0.0"),
                                          param("double", "yy", default_value = "0.0"),
                                          param("double", "yz", default_value = "0.0"),
                                          param("double", "zx", default_value = "0.0"),
                                          param("double", "zy", default_value = "0.0"),
                                          param("double", "zz", default_value = "0.0")])
        self.addSymTensorMethods(self.SymTensor1d, "SymTensor1d", 1)
        self.addSymTensorMethods(self.SymTensor2d, "SymTensor2d", 2)
        self.addSymTensorMethods(self.SymTensor3d, "SymTensor3d", 3)

        # Add methods to ThirdRankTensors.
        self.addThirdRankTensorMethods(self.ThirdRankTensor1d, 1)
        self.addThirdRankTensorMethods(self.ThirdRankTensor2d, 2)
        self.addThirdRankTensorMethods(self.ThirdRankTensor3d, 3)

        # Add methods to FourthRankTensors.
        self.addFourthRankTensorMethods(self.FourthRankTensor1d, 1)
        self.addFourthRankTensorMethods(self.FourthRankTensor2d, 2)
        self.addFourthRankTensorMethods(self.FourthRankTensor3d, 3)

        # Add methods to FifthRankTensors.
        self.addFifthRankTensorMethods(self.FifthRankTensor1d, 1)
        self.addFifthRankTensorMethods(self.FifthRankTensor2d, 2)
        self.addFifthRankTensorMethods(self.FifthRankTensor3d, 3)

        # Add methods to EigenStructs.
        self.addEigenStructMethods(self.EigenStruct1d, 1)
        self.addEigenStructMethods(self.EigenStruct2d, 2)
        self.addEigenStructMethods(self.EigenStruct3d, 3)

        # Add methods to planes.
        self.addPlaneMethods(self.Plane1d, 1)
        self.addPlaneMethods(self.Plane2d, 2)
        self.addPlaneMethods(self.Plane3d, 3)

        # Add methods to Box1d.
        self.addBox1dMethods(self.Box1d)

        # Add Polygon methods.
        self.addFacet2dMethods(self.Facet2d)
        self.addPolygonMethods(self.Polygon)
        self.space.add_function("aggregateFacetedVolumes", "Spheral::Polygon",
                                [constrefparam("vector_of_FacetedVolume2d", "loops")],
                                template_parameters = ["Spheral::Dim<2>"],
                                custom_name = "aggregateFacetedVolumes")

        # Add Polyhedron methods.
        self.addFacet3dMethods(self.Facet3d)
        self.addPolyhedronMethods(self.Polyhedron)
        self.space.add_function("aggregateFacetedVolumes", "Spheral::Polyhedron",
                                [constrefparam("vector_of_FacetedVolume3d", "loops")],
                                template_parameters = ["Spheral::Dim<3>"],
                                custom_name = "aggregateFacetedVolumes")

        generateStdVectorBindings(self.vector_of_Facet2d, "Spheral::Facet2d", "vector_of_Facet2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_Facet3d, "Spheral::Facet3d", "vector_of_Facet3d", indexAsPointer=True)

        generateStdVectorBindings(self.vector_of_FacetedVolume1d, "Spheral::Box1d", "vector_of_FacetedVolume1d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_FacetedVolume2d, "Spheral::Polygon", "vector_of_FacetedVolume2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_FacetedVolume3d, "Spheral::Polyhedron", "vector_of_FacetedVolume3d", indexAsPointer=True)

        generateStdVectorBindings(self.vector_of_vector_of_FacetedVolume1d, "vector_of_FacetedVolume1d", "vector_of_vector_of_FacetedVolume1d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_vector_of_FacetedVolume2d, "vector_of_FacetedVolume2d", "vector_of_vector_of_FacetedVolume2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_vector_of_FacetedVolume3d, "vector_of_FacetedVolume3d", "vector_of_vector_of_FacetedVolume3d", indexAsPointer=True)

        generateStdVectorBindings(self.vector_of_Plane1d, "Spheral::Plane1d", "vector_of_Plane1d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_Plane2d, "Spheral::Plane2d", "vector_of_Plane2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_Plane3d, "Spheral::Plane3d", "vector_of_Plane3d", indexAsPointer=True)

        # Add the free functions.
        for dim in self.dims:
            self.addDimFunctions(self.space, dim)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["Spheral"]

    #-------------------------------------------------------------------------------
    # Helper method for wrapping Vector.
    #-------------------------------------------------------------------------------
    def addVectorMethods(self, x, ndim):
    
        vecName = "Vector%id" % ndim
        vec3Name = "Vector3d"
        tenName = "Tensor%id" % ndim
        symTenName= "SymTensor%id" % ndim
        tensor_class = findObject(self.space, tenName)
    
        me = vecName

        # Instance attributes.
        x.add_static_attribute("nDimensions", "unsigned int", True)
        x.add_static_attribute("numElements", "unsigned int", True)
        x.add_static_attribute("zero", me, True)
        x.add_static_attribute("one", me, True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])
        x.add_constructor([param("double", "x")])
        x.add_constructor([param("double", "x"), param("double", "y")])
        x.add_constructor([param("double", "x"), param("double", "y"), param("double", "z")])
        x.add_constructor([param(me, "rhs")])
        x.add_function_as_constructor("constructGeomTypeFromSequence<%s>" % vecName,
                                      retval(ptr(vecName), caller_owns_return=True),
                                      [param("PyObject*", "sequence", transfer_ownership=False)])
    
        # x, y, z
        x.add_instance_attribute("x", "double", False, "x", "x")
        x.add_instance_attribute("y", "double", False, "y", "y")
        x.add_instance_attribute("z", "double", False, "z", "z")
    
        # Index by indices.
        x.add_method("operator()", "double", [param("int", "index")], custom_name="__call__")
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "unsigned int",
                                 [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__len__")
        x.add_function_as_method("indexGeomType", "double",
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [me],
                                 custom_name = "__getitem__")
        x.add_function_as_method("assignToGeomType", "int",
                                 [param(me, "self"), param("int", "index"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__setitem__")
        x.add_function_as_method("sliceGeomType", "vector_of_double",
                                 [param(me, "self"), param("int", "i1"), param("int", "i2")],
                                 template_parameters = [me],
                                 custom_name = "__getslice__")
        x.add_function_as_method("containsGeomType", "int",
                                 [param(me, "self"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__contains__")

        x.add_method("Zero", None, [])
    
        # Operators.
        x.add_unary_numeric_operator("-")
    
        x.add_binary_numeric_operator("+")
        x.add_binary_numeric_operator("-")
        x.add_binary_numeric_operator("*", result_cppclass = tensor_class)
    
        x.add_inplace_numeric_operator("+=")
        x.add_inplace_numeric_operator("-=")
    
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        x.add_binary_numeric_operator("*", left_cppclass = double)
    
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        # Comparisons.
        x.add_method("compare", "int", [param(me, "other")], is_const=True)
        x.add_method("compare", "int", [param("double", "other")], is_const=True)
    
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")
    
        # Misc methods.
        x.add_method("dot", "double", [param(me, "other")], is_const=True)
        x.add_method("cross", vec3Name, [param(me, "other")], is_const=True)
        x.add_method("dyad", tenName, [param(me, "other")], is_const=True)
        x.add_method("selfdyad", symTenName, [], is_const=True)
        x.add_method("unitVector", me, [], is_const=True)
        x.add_method("magnitude", "double", [], is_const=True)
        x.add_method("magnitude2", "double", [], is_const=True)
        x.add_method("minElement", "double", [], is_const=True)
        x.add_method("maxElement", "double", [], is_const=True)
        x.add_method("maxAbsElement", "double", [], is_const=True)
        x.add_method("sumElements", "double", [], is_const=True)
    
        # String representations.
        x.add_function_as_method("printReprVector", "std::string", [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__repr__")
        x.add_output_stream_operator()
    
        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping Geom3Vector.
    #-------------------------------------------------------------------------------
    def addGeom3VectorMethods(self, x):
    
        me = "Geom3Vector"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping Tensor.
    #-------------------------------------------------------------------------------
    def addTensorMethods(self, x, me, ndim):
    
        vec = "Vector%id" % ndim
        ten = "Tensor%id" % ndim
        symten = "SymTensor%id" % ndim

        vecType = findObject(self.space, vec)
        tenType = findObject(self.space, ten)
        symtenType = findObject(self.space, symten)
    
        # Instance attributes.
        x.add_static_attribute("nDimensions", "unsigned int", True)
        x.add_static_attribute("numElements", "unsigned int", True)
        x.add_static_attribute("zero", me, True)
        x.add_static_attribute("one", me, True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])
        x.add_constructor([param(ten, "rhs")])
        x.add_constructor([param(symten, "rhs")])
        x.add_function_as_constructor("constructGeomTypeFromSequence<%s>" % me,
                                      retval(ptr(me), caller_owns_return=True),
                                      [param("PyObject*", "sequence", transfer_ownership=False)])

        # Components.
        x.add_instance_attribute("xx", "double", False, "xx", "xx")
        x.add_instance_attribute("xy", "double", False, "xy", "xy")
        x.add_instance_attribute("xz", "double", False, "xz", "xz")
        x.add_instance_attribute("yx", "double", False, "yx", "yx")
        x.add_instance_attribute("yy", "double", False, "yy", "yy")
        x.add_instance_attribute("yz", "double", False, "yz", "yz")
        x.add_instance_attribute("zx", "double", False, "zx", "zx")
        x.add_instance_attribute("zy", "double", False, "zy", "zy")
        x.add_instance_attribute("zz", "double", False, "zz", "zz")
    
        # Index by indices.
        x.add_method("operator()", "double", [param("int", "row"), param("int", "column")], custom_name="__call__")
        x.add_function_as_method("assignSecondRankTensorElement",
                                 None,
                                 [param(me, "self"), param("int", "i"), param("int", "j"), param("double", "val")],
                                 template_parameters = [me],
                                 custom_name = "__call__")
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "unsigned int",
                                 [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__len__")
        x.add_function_as_method("indexGeomType", "double",
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [me],
                                 custom_name = "__getitem__")
        x.add_function_as_method("assignToGeomType", "int",
                                 [param(me, "self"), param("int", "index"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__setitem__")
        x.add_function_as_method("sliceGeomType", "vector_of_double",
                                 [param(me, "self"), param("int", "i1"), param("int", "i2")],
                                 template_parameters = [me],
                                 custom_name = "__getslice__")
        x.add_function_as_method("containsGeomType", "int",
                                 [param(me, "self"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__contains__")

        # Extract rows/cols.
        x.add_method("getRow", vec, [param("int", "row")], is_const=True)
        x.add_method("getColumn", vec, [param("int", "column")], is_const=True)
    
        x.add_method("Zero", None, [])
        x.add_method("Identity", None, [])
    
        # Operators.
        x.add_unary_numeric_operator("-")
    
        x.add_inplace_numeric_operator("+=")
        x.add_inplace_numeric_operator("-=")
    
        #x.add_inplace_numeric_operator("+=", right = "double")
        #x.add_inplace_numeric_operator("-=", right = "double")
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        x.add_binary_numeric_operator("+", result_cppclass = tenType, right = ten)
        x.add_binary_numeric_operator("-", result_cppclass = tenType, right = ten)
        x.add_binary_numeric_operator("*", result_cppclass = tenType, right = ten)
    
        x.add_binary_numeric_operator("+", right = symten)
        x.add_binary_numeric_operator("-", right = symten)
        x.add_binary_numeric_operator("*", result_cppclass = tenType, right = symten)
        x.add_binary_numeric_operator("*", result_cppclass = vecType, right = vec)
    
        #x.add_binary_numeric_operator("+", right = "double")
        #x.add_binary_numeric_operator("-", right = "double")
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        x.add_binary_numeric_operator("*", left_cppclass = double)
    
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")
    
        # Misc methods.
        x.add_method("Symmetric", symten, [], is_const=True)
        x.add_method("SkewSymmetric", ten, [], is_const=True)
        x.add_method("Transpose", me, [], is_const=True)
        x.add_method("Inverse", me, [], is_const=True)
        x.add_method("diagonalElements", vec, [], is_const=True)
        x.add_method("Trace", "double", [], is_const=True)
        x.add_method("Determinant", "double", [], is_const=True)
        x.add_method("dot", vec, [param(vec, "other")], is_const=True)
        x.add_method("dot", ten, [param(ten, "other")], is_const=True)
        x.add_method("dot", ten, [param(symten, "other")], is_const=True)
        x.add_method("doubledot", "double", [param(ten, "other")], is_const=True)
        x.add_method("doubledot", "double", [param(symten, "other")], is_const=True)
        x.add_method("selfDoubledot", "double", [], is_const=True)
        x.add_method("square", me, [], is_const=True)
        x.add_method("squareElements", me, [], is_const=True)
        x.add_method("eigenValues", vec, [], is_const=True)
        x.add_method("rotationalTransform", None, [param(ten, "R")])
        x.add_method("maxAbsElement", "double", [], is_const=True)
    
        x.add_function_as_method("printReprTensor", "std::string", [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__repr__")
        x.add_output_stream_operator()
    
        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping SymTensor.
    #-------------------------------------------------------------------------------
    def addSymTensorMethods(self, x, me, ndim):
    
        vec = "Vector%id" % ndim
        ten = "Tensor%id" % ndim
        symten = "SymTensor%id" % ndim
        eigenstruct = "EigenStruct%id" % ndim

        vecType = findObject(self.space, vec)
        tenType = findObject(self.space, ten)
        symtenType = findObject(self.space, symten)
    
        # First add the generic Tensor methods.
        self.addTensorMethods(x, me, ndim)

        # The few unique methods to SymTensors.
        x.add_method("cube", me, [], is_const=True)
        x.add_method("sqrt", me, [], is_const=True)
        x.add_method("cuberoot", me, [], is_const=True)
        x.add_method("pow", me, [param("double", "p")], is_const=True)
        x.add_method("eigenVectors", eigenstruct, [], is_const=True)
    
        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping ThirdRankTensor.
    #-------------------------------------------------------------------------------
    def addThirdRankTensorMethods(self, x, ndim):
    
        me = "ThirdRankTensor%id" % ndim

        # Add the base methods.
        self.addRankNTensorMethods(x, ndim, me)

        # Index by indices.
        x.add_method("operator()", "double", [param("int", "i"), param("int", "j"), param("int", "k")], custom_name="__call__")
        x.add_function_as_method("assignThirdRankTensorElement",
                                 None,
                                 [param(me, "self"), param("int", "i"), param("int", "j"), param("int", "k"), param("double", "val")],
                                 template_parameters = [me],
                                 custom_name = "__call__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping FourthRankTensor.
    #-------------------------------------------------------------------------------
    def addFourthRankTensorMethods(self, x, ndim):
    
        me = "FourthRankTensor%id" % ndim

        # Add the base methods.
        self.addRankNTensorMethods(x, ndim, me)

        # Index by indices.
        x.add_method("operator()", "double", [param("int", "i"), param("int", "j"), param("int", "k"), param("int", "m")], custom_name="__call__")
        x.add_function_as_method("assignFourthRankTensorElement",
                                 None,
                                 [param(me, "self"), param("int", "i"), param("int", "j"), param("int", "k"), param("int", "m"), param("double", "val")],
                                 template_parameters = [me],
                                 custom_name = "__call__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping FifthRankTensor.
    #-------------------------------------------------------------------------------
    def addFifthRankTensorMethods(self, x, ndim):
    
        me = "FifthRankTensor%id" % ndim

        # Add the base methods.
        self.addRankNTensorMethods(x, ndim, me)

        # Index by indices.
        x.add_method("operator()", "double", [param("int", "i"), param("int", "j"), param("int", "k"), param("int", "m"), param("int", "n")], custom_name="__call__")
        x.add_function_as_method("assignFifthRankTensorElement",
                                 None,
                                 [param(me, "self"), param("int", "i"), param("int", "j"), param("int", "k"), param("int", "m"), param("int", "n"), param("double", "val")],
                                 template_parameters = [me],
                                 custom_name = "__call__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping methods common to  RankNTensor descendants.
    #-------------------------------------------------------------------------------
    def addRankNTensorMethods(self, x, ndim, me):
    
        # Instance attributes.
        x.add_static_attribute("nrank", "unsigned int", is_const=True)
        x.add_static_attribute("nDimensions", "unsigned int", is_const=True)
        x.add_static_attribute("numElements", "unsigned int", is_const=True)
        x.add_static_attribute("zero", me, True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(me, "rhs")])
        x.add_constructor([param("double", "val")])
        x.add_constructor([param(me, "rhs")])
        x.add_function_as_constructor("constructGeomTypeFromSequence<%s>" % me,
                                      retval(ptr(me), caller_owns_return=True),
                                      [param("PyObject*", "sequence", transfer_ownership=False)])
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "unsigned int",
                                 [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__len__")
        x.add_function_as_method("indexGeomType", "double",
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [me],
                                 custom_name = "__getitem__")
        x.add_function_as_method("assignToGeomType", "int",
                                 [param(me, "self"), param("int", "index"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__setitem__")
        x.add_function_as_method("sliceGeomType", "vector_of_double",
                                 [param(me, "self"), param("int", "i1"), param("int", "i2")],
                                 template_parameters = [me],
                                 custom_name = "__getslice__")
        x.add_function_as_method("containsGeomType", "int",
                                 [param(me, "self"), param("double", "value")],
                                 template_parameters = [me],
                                 custom_name = "__contains__")
    
        x.add_method("Zero", None, [])
    
        # Operators.
        x.add_unary_numeric_operator("-")
    
        x.add_inplace_numeric_operator("+=")
        x.add_inplace_numeric_operator("-=")
    
        x.add_binary_numeric_operator("+")
        x.add_binary_numeric_operator("-")
    
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")
    
        # Misc methods.
        x.add_method("doubledot", "double", [param(me, "rhs")], is_const=True)
        x.add_method("squareElements", me, [], is_const=True)
        x.add_method("maxAbsElement", "double", [], is_const=True)
        x.add_output_stream_operator()

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping EigenStruct.
    #-------------------------------------------------------------------------------
    def addEigenStructMethods(self, x, ndim):
    
        me = "EigenStruct%id" % ndim
        vec = "Vector%id" % ndim
        ten = "Tensor%id" % ndim

        # Instance attributes.
        x.add_instance_attribute("eigenValues", vec)
        x.add_instance_attribute("eigenVectors", ten)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([param(me, "rhs")])
    
        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping GeomPlane.
    #-------------------------------------------------------------------------------
    def addPlaneMethods(self, x, ndim):

        me = "Plane%id" % ndim
        vec = "Vector%id" % ndim
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([param(me, "rhs")])
        x.add_constructor([param(vec, "point"), param(vec, "normal")])
        x.add_constructor([param("vector_of_%s" % vec, "points")])
    
        # Attributes.
        x.add_instance_attribute("point", vec, is_const = False, getter = "point", setter = "point")
        x.add_instance_attribute("normal", vec, is_const = False, getter = "normal", setter = "normal")
    
        # Operators.
        x.add_unary_numeric_operator("-")

        # Methods.
        x.add_method("signedDistance", "double", [param(vec, "point")], is_const = True)
        x.add_method("minimumDistance", "double", [param(vec, "point")], is_const = True)
        x.add_method("closestPointOnPlane", vec, [param(vec, "point")], is_const = True)
        x.add_method("parallel", "bool", [param(me, "rhs")], is_const = True)
        x.add_method("valid", "bool", [], is_const = True)
        x.add_method("compare", "int", [constrefparam(vec, "point")], is_const=True)
        
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
    
        return

    #-------------------------------------------------------------------------------
    # Helper method for Box1d.
    #-------------------------------------------------------------------------------
    def addBox1dMethods(self, x):

        me = "Spheral::Box1d"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam("vector_of_Vector1d", "points")])
        x.add_constructor([param("Vector1d", "center"), param("double", "extent")])
        x.add_constructor([constrefparam(me, "rhs")])
    
        # Attributes.
        x.add_instance_attribute("center", "Vector1d", getter = "center", setter = "center")
        x.add_instance_attribute("extent", "double", getter = "extent", setter = "extent")
        x.add_instance_attribute("xmin", "Vector1d", getter = "xmin", is_const = True)
        x.add_instance_attribute("xmax", "Vector1d", getter = "xmax", is_const = True)
        x.add_instance_attribute("volume", "double", getter="volume", is_const=True)
    
        # Methods.
        x.add_method("contains", "bool", [constrefparam("Vector1d", "point"),
                                          param("bool", "countBoundary", default_value="true"),
                                          param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("convexContains", "bool", [constrefparam("Vector1d", "point"),
                                                param("bool", "countBoundary", default_value="true"),
                                                param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("convexIntersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("distance", "double", [constrefparam("Vector1d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector1d", [constrefparam("Vector1d", "point")], is_const=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector1d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector1d", "&%s::vertices" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertices")
    
        x.add_inplace_numeric_operator("+=", right = "Vector1d")
        x.add_inplace_numeric_operator("-=", right = "Vector1d")
        x.add_binary_numeric_operator("+", right = "Vector1d")
        x.add_binary_numeric_operator("-", right = "Vector1d")

        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        return

    #-------------------------------------------------------------------------------
    # Facet2d
    #-------------------------------------------------------------------------------
    def addFacet2dMethods(self, x):

        me = "Spheral::Facet2d"

        # Constructors.
        x.add_constructor([constrefparam("vector_of_Vector2d", "vertices"),
                           param("unsigned int", "point1"),
                           param("unsigned int", "point2")])
        x.add_constructor([refparam(me, "rhs")])
    
        # Methods.
        x.add_method("compare", "int", [constrefparam("Vector2d", "point"),
                                        param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("compare", "int", [constrefparam("vector_of_Vector2d", "points"),
                                        param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("distance", "double", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector2d", [constrefparam("Vector2d", "point")], is_const=True)

        # Attributes.
        x.add_instance_attribute("point1", "Vector2d", getter="point1", is_const=True)
        x.add_instance_attribute("point2", "Vector2d", getter="point2", is_const=True)
        x.add_instance_attribute("ipoint1", "unsigned int", getter="ipoint1", is_const=True)
        x.add_instance_attribute("ipoint2", "unsigned int", getter="ipoint2", is_const=True)
        x.add_instance_attribute("ipoints", "vector_of_unsigned", getter="ipoints", is_const=True)
        x.add_instance_attribute("normal", "Vector2d", getter="normal", is_const=True)
        x.add_instance_attribute("position", "Vector2d", getter="position", is_const=True)
        x.add_instance_attribute("area", "double", getter="area", is_const=True)

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

    #-------------------------------------------------------------------------------
    # Facet3d
    #-------------------------------------------------------------------------------
    def addFacet3dMethods(self, x):

        me = "Spheral::Facet3d"

        # Constructors.
        x.add_constructor([constrefparam("vector_of_Vector3d", "vertices"),
                           constrefparam("vector_of_unsigned", "ipoints"),
                           constrefparam("Vector3d", "normal")])
        x.add_constructor([refparam(me, "rhs")])

        # Methods.
        x.add_method("compare", "int", [constrefparam("Vector3d", "point"),
                                        param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("compare", "int", [constrefparam("vector_of_Vector3d", "point"),
                                        param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("point", "Vector3d", [param("unsigned int", "index")], is_const=True)
        x.add_method("distance", "double", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector3d", [constrefparam("Vector3d", "point")], is_const=True)

        # Attributes.
        x.add_instance_attribute("ipoints", "vector_of_unsigned", getter="ipoints", is_const=True)
        x.add_instance_attribute("normal", "Vector3d", getter="normal", is_const=True)
        x.add_instance_attribute("position", "Vector3d", getter="position", is_const=True)
        x.add_instance_attribute("area", "double", getter="area", is_const=True)

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

    #-------------------------------------------------------------------------------
    # Polygon
    #-------------------------------------------------------------------------------
    def addPolygonMethods(self, x):

        me = "Spheral::Polygon"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam("vector_of_Vector2d", "points")])
        x.add_constructor([constrefparam("vector_of_Vector2d", "points"),
                           constrefparam("vector_of_vector_of_unsigned", "facets")])
        x.add_constructor([constrefparam(me, "rhs")])
    
        # Methods.
        x.add_method("contains", "bool", [constrefparam("Vector2d", "point"),
                                          param("bool", "countBoundary", default_value="true"),
                                          param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("convexContains", "bool", [constrefparam("Vector2d", "point"),
                                                param("bool", "countBoundary", default_value="true"),
                                                param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("convexIntersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam("pair_Vector2d_Vector2d", "rhs")], is_const=True)
        x.add_method("intersect", "vector_of_Vector2d", [constrefparam("Vector2d", "x0"),
                                                         constrefparam("Vector2d", "x1")], is_const=True)
        x.add_method("centroid", "Vector2d", [], is_const=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector2d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector2d", "&%s::vertices" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertices")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Facet2d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Facet2d", "&%s::facets" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "facets")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector2d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector2d", "&%s::vertexUnitNorms" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertexUnitNorms")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_vector_of_unsigned"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_vector_of_unsigned", "&%s::vertexFacetConnectivity" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertexFacetConnectivity")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_vector_of_unsigned"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_vector_of_unsigned", "&%s::facetFacetConnectivity" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "facetFacetConnectivity")
        x.add_method("closestFacet", "unsigned int", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("distance", "double", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector2d", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("reconstruct", None, [constrefparam("vector_of_Vector2d", "vertices"),
                                           constrefparam("vector_of_vector_of_unsigned", "facetVertices")])
        x.add_method("convex", "bool", [param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("setBoundingBox", None, [])

        x.add_inplace_numeric_operator("+=", right = "Vector2d")
        x.add_inplace_numeric_operator("-=", right = "Vector2d")
        x.add_binary_numeric_operator("+", right = "Vector2d")
        x.add_binary_numeric_operator("-", right = "Vector2d")

        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")

        # Attributes.
        x.add_instance_attribute("xmin", "Vector2d", getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", "Vector2d", getter="xmax", is_const=True)
        x.add_instance_attribute("edges", "vector_of_pair_unsigned_unsigned", getter="edges", is_const=True)
        x.add_instance_attribute("facetVertices", "vector_of_vector_of_unsigned", getter="facetVertices", is_const=True)
        x.add_instance_attribute("volume", "double", getter="volume", is_const=True)

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        # String representation.
        x.add_output_stream_operator()

        return

    #-------------------------------------------------------------------------------
    # Polyhedron
    #-------------------------------------------------------------------------------
    def addPolyhedronMethods(self, x):

        me = "Spheral::Polyhedron"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam("vector_of_Vector3d", "points")])
        x.add_constructor([refparam("vector_of_Vector3d", "points"),
                           refparam("vector_of_vector_of_unsigned", "facetIndices")])
        x.add_constructor([constrefparam(me, "rhs")])
    
        # Methods.
        x.add_method("contains", "bool", [constrefparam("Vector3d", "point"),
                                          param("bool", "countBoundary", default_value="true"),
                                          param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("convexContains", "bool", [constrefparam("Vector3d", "point"),
                                                param("bool", "countBoundary", default_value="true"),
                                                param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("convexIntersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam("pair_Vector3d_Vector3d", "rhs")], is_const=True)
        x.add_method("centroid", "Vector3d", [], is_const=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector3d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector3d", "&%s::vertices" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertices")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Facet3d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Facet3d", "&%s::facets" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "facets")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector3d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector3d", "&%s::vertexUnitNorms" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertexUnitNorms")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_vector_of_unsigned"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_vector_of_unsigned", "&%s::vertexFacetConnectivity" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertexFacetConnectivity")
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_vector_of_unsigned"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_vector_of_unsigned", "&%s::facetFacetConnectivity" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "facetFacetConnectivity")
        x.add_method("closestFacet", "unsigned int", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("distance", "double", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector3d", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("reconstruct", None, [constrefparam("vector_of_Vector3d", "vertices"),
                                           constrefparam("vector_of_vector_of_unsigned", "facetVertices"),
                                           constrefparam("vector_of_Vector3d", "facetNormals")])
        x.add_method("convex", "bool", [param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("setBoundingBox", None, [])

        x.add_inplace_numeric_operator("+=", right = "Vector3d")
        x.add_inplace_numeric_operator("-=", right = "Vector3d")
        x.add_binary_numeric_operator("+", right = "Vector3d")
        x.add_binary_numeric_operator("-", right = "Vector3d")

        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")

        # Attributes.
        x.add_instance_attribute("xmin", "Vector3d", getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", "Vector3d", getter="xmax", is_const=True)
        x.add_instance_attribute("edges", "vector_of_pair_unsigned_unsigned", getter="edges", is_const=True)
        x.add_instance_attribute("facetVertices", "vector_of_vector_of_unsigned", getter="facetVertices", is_const=True)
        x.add_instance_attribute("facetNormals", "vector_of_Vector3d", getter="facetNormals", is_const=True)
        x.add_instance_attribute("volume", "double", getter="volume", is_const=True)

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        # String representation.
        x.add_output_stream_operator()

        return

    #-------------------------------------------------------------------------------
    # Add the free functions for a given dimensionality.
    #-------------------------------------------------------------------------------
    def addDimFunctions(self, space, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        thirdranktensor = "ThirdRankTensor%id" % ndim
        fourthranktensor = "FourthRankTensor%id" % ndim
        fifthranktensor = "FifthRankTensor%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim

        space.add_function("computeEigenValues", None, 
                           [constrefparam(symtensorfield, "field"),
                            refparam(vectorfield, "eigenValues"),
                            refparam(tensorfield, "eigenVectors")],
                           template_parameters = ["Dim<%i> " % ndim],
                           custom_name = "computeEigenValues",
                           docstring = "Compute the eigen values & vectors for a field of symmetric tensors.")

        space.add_function("invertRankNTensor", fourthranktensor,
                           [constrefparam(fourthranktensor, "tensor")],
                           template_parameters = [fourthranktensor],
                           custom_name = "invertRankNTensor",
                           docstring = "Invert a fouth rank tensor.")

        # Inner & outer products with doubles
        for type in (vector, tensor, symtensor, thirdranktensor, fourthranktensor, fifthranktensor):
            space.add_function("innerProduct", type, 
                               [param("double", "lhs"), constrefparam(type, "rhs")], 
                               foreign_cpp_namespace = "Spheral",
                               template_parameters=[type],
                               custom_name="innerProduct")
            space.add_function("innerProduct", type,
                               [constrefparam(type, "lhs"), param("double", "rhs")], 
                               foreign_cpp_namespace = "Spheral",
                               template_parameters=[type], 
                               custom_name="innerProduct")
            space.add_function("outerProduct", type, 
                               [param("double", "lhs"), constrefparam(type, "rhs")], 
                               foreign_cpp_namespace = "Spheral",
                               template_parameters=[type], 
                               custom_name="outerProduct")
            space.add_function("outerProduct", type, 
                               [constrefparam(type, "lhs"), param("double", "rhs")], 
                               foreign_cpp_namespace = "Spheral",
                               template_parameters=[type],
                               custom_name="outerProduct")

        # vec . vec
        space.add_function("innerProduct", "double", 
                           [constrefparam(vector, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # tensor . vec
        space.add_function("innerProduct", vector, 
                           [constrefparam(tensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", vector, 
                           [constrefparam(symtensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # vec . tensor
        space.add_function("innerProduct", vector, 
                           [constrefparam(vector, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", vector, 
                           [constrefparam(vector, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # thirdranktensor . vector
        space.add_function("innerProduct", tensor, 
                           [constrefparam(thirdranktensor, "rhs"), constrefparam(vector, "lhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", tensor, 
                           [constrefparam(vector, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # tensor . tensor
        space.add_function("innerProduct", tensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", tensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", tensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", tensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # tensor . thirdranktensor
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # thirdranktensor . thirdranktensor
        space.add_function("innerProduct", fourthranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # fourthranktensor . vector
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", thirdranktensor, 
                           [constrefparam(vector, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # fourthranktensor . tensor
        space.add_function("innerProduct", fourthranktensor, 
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", fourthranktensor, 
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", fourthranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", fourthranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # fourthranktensor . thirdranktensor
        space.add_function("innerProduct", fifthranktensor, 
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")
        space.add_function("innerProduct", fifthranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerProduct")

        # vec outer vec
        space.add_function("outerProduct", tensor, 
                           [constrefparam(vector, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # tensor outer vec
        space.add_function("outerProduct", thirdranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", thirdranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", thirdranktensor, 
                           [constrefparam(vector, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", thirdranktensor, 
                           [constrefparam(vector, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # thirdranktensor outer vec
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(vector, "rhs"), constrefparam(thirdranktensor, "lhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # tensor outer tensor
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fourthranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # fourthranktensor outer vec
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(vector, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(vector, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # thirdranktensor outer tensor
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(tensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")
        space.add_function("outerProduct", fifthranktensor, 
                           [constrefparam(symtensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="outerProduct")

        # tensor .. tensor
        space.add_function("innerDoubleProduct", "double", 
                           [constrefparam(tensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", "double", 
                           [constrefparam(tensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", "double", 
                           [constrefparam(symtensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", "double", 
                           [constrefparam(symtensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # tensor .. thirdranktensor
        space.add_function("innerDoubleProduct", vector,
                           [constrefparam(tensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", vector,
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", vector,
                           [constrefparam(symtensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", vector,
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # thirdranktensor .. thirdranktensor
        space.add_function("innerDoubleProduct", tensor,
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # tensor .. fourthranktensor
        space.add_function("innerDoubleProduct", tensor,
                           [constrefparam(tensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", tensor,
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", tensor,
                           [constrefparam(symtensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", tensor,
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # thirdranktensor .. fourthranktensor
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # fourthranktensor .. fourthranktensor
        space.add_function("innerDoubleProduct", fourthranktensor,
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # tensor .. fifthranktensor
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(tensor, "lhs"), constrefparam(fifthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(fifthranktensor, "lhs"), constrefparam(tensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(symtensor, "lhs"), constrefparam(fifthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", thirdranktensor,
                           [constrefparam(fifthranktensor, "lhs"), constrefparam(symtensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # thirdranktensor .. fifthranktensor
        space.add_function("innerDoubleProduct", fourthranktensor,
                           [constrefparam(thirdranktensor, "lhs"), constrefparam(fifthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", fourthranktensor,
                           [constrefparam(fifthranktensor, "lhs"), constrefparam(thirdranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        # fourthranktensor .. fifthranktensor
        space.add_function("innerDoubleProduct", fifthranktensor,
                           [constrefparam(fourthranktensor, "lhs"), constrefparam(fifthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")
        space.add_function("innerDoubleProduct", fifthranktensor,
                           [constrefparam(fifthranktensor, "lhs"), constrefparam(fourthranktensor, "rhs")],
                           foreign_cpp_namespace = "Spheral",
                           template_parameters=[dim],
                           custom_name="innerDoubleProduct")

        return
