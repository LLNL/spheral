from pybindgen import *

import sys
sys.path.append("CXXTypes")
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
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Geometry/GeometryTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral

        # Expose types.
        self.Vector1d = addObject(Spheral, "Vector1d")
        self.Vector2d = addObject(Spheral, "Vector2d")
        self.Vector3d = addObject(Spheral, "Vector3d")

        self.Geom3Vector = addObject(Spheral, "Geom3Vector")
        
        self.Tensor1d = addObject(Spheral, "Tensor1d")
        self.Tensor2d = addObject(Spheral, "Tensor2d")
        self.Tensor3d = addObject(Spheral, "Tensor3d")

        self.SymTensor1d = addObject(Spheral, "SymTensor1d")
        self.SymTensor2d = addObject(Spheral, "SymTensor2d")
        self.SymTensor3d = addObject(Spheral, "SymTensor3d")

        self.ThirdRankTensor1d = addObject(Spheral, "ThirdRankTensor1d")
        self.ThirdRankTensor2d = addObject(Spheral, "ThirdRankTensor2d")
        self.ThirdRankTensor3d = addObject(Spheral, "ThirdRankTensor3d")

        self.EigenStruct1d = addObject(Spheral, "EigenStruct1d")
        self.EigenStruct2d = addObject(Spheral, "EigenStruct2d")
        self.EigenStruct3d = addObject(Spheral, "EigenStruct3d")

        self.Plane1d = addObject(Spheral, "Plane1d")
        self.Plane2d = addObject(Spheral, "Plane2d")
        self.Plane3d = addObject(Spheral, "Plane3d")

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

        # Add methods to Tensors.
        self.addTensorMethods(self.Tensor1d, "Tensor1d", 1)
        self.addTensorMethods(self.Tensor2d, "Tensor2d", 2)
        self.addTensorMethods(self.Tensor3d, "Tensor3d", 3)
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

        # Add methods to SymTensors.
        self.addSymTensorMethods(self.SymTensor1d, "SymTensor1d", 1)
        self.addSymTensorMethods(self.SymTensor2d, "SymTensor2d", 2)
        self.addSymTensorMethods(self.SymTensor3d, "SymTensor3d", 3)
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

        # Add methods to ThirdRankTensors.
        self.addThirdRankTensorMethods(self.ThirdRankTensor1d, 1)
        self.addThirdRankTensorMethods(self.ThirdRankTensor2d, 2)
        self.addThirdRankTensorMethods(self.ThirdRankTensor3d, 3)

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

        # Add Polyhedron methods.
        self.addFacet3dMethods(self.Facet3d)
        self.addPolyhedronMethods(self.Polyhedron)

        generateStdVectorBindings(self.vector_of_Facet2d, "Spheral::Facet2d", "vector_of_Facet2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_Facet3d, "Spheral::Facet3d", "vector_of_Facet3d", indexAsPointer=True)

        # Add the free functions.
        self.addDimFunctions(self.space, 1)
        self.addDimFunctions(self.space, 2)
        self.addDimFunctions(self.space, 3)

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
        tensor_class = self.space.wrapObjs[tenName]
    
        me = vecName

        # Instance attributes.
        x.add_static_attribute("nDimensions", "int", True)
        x.add_static_attribute("zero", me, True)
        x.add_static_attribute("one", me, True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("double", "x")])
        x.add_constructor([param("double", "x"), param("double", "y")])
        x.add_constructor([param("double", "x"), param("double", "y"), param("double", "z")])
        x.add_constructor([param(me, "rhs")])
    
        # x, y, z
        x.add_instance_attribute("x", "double", False, "x", "x")
        x.add_instance_attribute("y", "double", False, "y", "y")
        x.add_instance_attribute("z", "double", False, "z", "z")
    
        # Index by indicies.
        x.add_method("operator()", "double", [param("int", "index")], custom_name="__call__")
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "int",
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
    
        x.add_binary_numeric_operator("+", right = "double")
        x.add_binary_numeric_operator("-", right = "double")
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        x.add_binary_numeric_operator("*", left_cppclass = double)
    
        x.add_inplace_numeric_operator("+=", right = "double")
        x.add_inplace_numeric_operator("-=", right = "double")
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
    # Helper method for wrapping Tensor.
    #-------------------------------------------------------------------------------
    def addTensorMethods(self, x, me, ndim):
    
        vec = "Vector%id" % ndim
        ten = "Tensor%id" % ndim
        symten = "SymTensor%id" % ndim

        vecType = self.space.wrapObjs[vec]
        tenType = self.space.wrapObjs[ten]
        symtenType = self.space.wrapObjs[symten]
    
        # Instance attributes.
        x.add_static_attribute("nDimensions", "int", True)
        x.add_static_attribute("zero", me, True)
        x.add_static_attribute("one", me, True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([param(ten, "rhs")])
        x.add_constructor([param(symten, "rhs")])
    
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
    
        # Index by indicies.
        x.add_method("operator()", "double", [param("int", "row"), param("int", "column")], custom_name="__call__")
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "int",
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
    
        x.add_inplace_numeric_operator("+=", right = "double")
        x.add_inplace_numeric_operator("-=", right = "double")
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        x.add_binary_numeric_operator("+", result_cppclass = tenType, right = ten)
        x.add_binary_numeric_operator("-", result_cppclass = tenType, right = ten)
        x.add_binary_numeric_operator("*", result_cppclass = tenType, right = ten)
    
        x.add_binary_numeric_operator("+", right = symten)
        x.add_binary_numeric_operator("-", right = symten)
        x.add_binary_numeric_operator("*", result_cppclass = tenType, right = symten)
        x.add_binary_numeric_operator("*", result_cppclass = vecType, right = vec)
    
        x.add_binary_numeric_operator("+", right = "double")
        x.add_binary_numeric_operator("-", right = "double")
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

        vecType = self.space.wrapObjs[vec]
        tenType = self.space.wrapObjs[ten]
        symtenType = self.space.wrapObjs[symten]
    
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

        # Instance attributes.
        x.add_static_attribute("nDimensions", "int", True)
        x.add_static_attribute("numElements", "int", True)
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("double", "val")])
        x.add_constructor([param(me, "rhs")])
    
        # Index by indicies.
        x.add_method("operator()", "double", [param("int", "i"), param("int", "j"), param("int", "k")], custom_name="__call__")
        x.add_function_as_method("assignThirdRankTensorElement",
                                 None,
                                 [param(me, "self"), param("int", "i"), param("int", "j"), param("int", "k"), param("double", "val")],
                                 template_parameters = [me],
                                 custom_name = "__call__")
    
        # Add sequence methods.
        x.add_function_as_method("sizeGeomType", "int",
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
    
        x.add_function_as_method("printReprThirdRankTensor", "std::string", [param(me, "self")],
                                 template_parameters = [me],
                                 custom_name = "__repr__")
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
    
        # Attributes.
        x.add_instance_attribute("point", vec, is_const = False, getter = "point", setter = "point")
        x.add_instance_attribute("normal", vec, is_const = False, getter = "normal", setter = "normal")
    
        # Operators.
        x.add_unary_numeric_operator("-")

        # Methods.
        x.add_method("minimumDistance", "double", [param(vec, "point")], is_const = True)
        x.add_method("parallel", "bool", [param(me, "rhs")], is_const = True)
        x.add_method("valid", "bool", [], is_const = True)
        x.add_method("compare", "int", [constrefparam(vec, "point")], is_const=True)
        
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
    
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
    
        # Methods.
        x.add_method("contains", "bool", [constrefparam("Vector1d", "point"),
                                          param("bool", "countBoundary", default_value="true"),
                                          param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("convexContains", "bool", [constrefparam("Vector1d", "point"),
                                                param("bool", "countBoundary", default_value="true"),
                                                param("double", "tol", default_value="1.0e-8")], is_const=True)
        x.add_method("intersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_method("convexIntersect", "bool", [constrefparam(me, "rhs")], is_const=True)
        x.add_function_as_method("const_reference_as_pointer",
                                 retval(ptr("vector_of_Vector1d"), reference_existing_object=True),
                                 [param(me, "self")],
                                 template_parameters = [me, "vector_of_Vector1d", "&%s::vertices" % me],
                                 foreign_cpp_namespace = "Spheral",
                                 custom_name = "vertices")
    
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
        box = "Wm5::WMBox2d"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam("vector_of_Vector2d", "points")])
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
        x.add_method("intersect", "bool", [constrefparam(box, "box")], is_const=True)
        x.add_method("intersect", "vector_of_Vector2d", [constrefparam("Vector2d", "x0"),
                                                         constrefparam("Vector2d", "x1")], is_const=True)
        x.add_method("centroid", "Vector2d", [], is_const=True)
#         x.add_method("vertices", retval("const vector_of_Vector2d&", reference_existing_object=True), [], is_const=True)
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
        x.add_method("distance", "double", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector2d", [constrefparam("Vector2d", "point")], is_const=True)
        x.add_method("reconstruct", None, [constrefparam("vector_of_Vector2d", "vertices"),
                                           constrefparam("vector_of_vector_of_unsigned", "facetVertices")])
        x.add_method("convex", "bool", [param("double", "tol", default_value="1.0e-8")], is_const=True)

        # Attributes.
        x.add_instance_attribute("xmin", "Vector2d", getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", "Vector2d", getter="xmax", is_const=True)
        x.add_instance_attribute("edges", "vector_of_pair_unsigned_unsigned", getter="edges", is_const=True)
        x.add_instance_attribute("facetVertices", "vector_of_vector_of_unsigned", getter="facetVertices", is_const=True)
        x.add_instance_attribute("volume", "double", getter="volume", is_const=True)

        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")

        return

    #-------------------------------------------------------------------------------
    # Polyhedron
    #-------------------------------------------------------------------------------
    def addPolyhedronMethods(self, x):

        me = "Spheral::Polyhedron"
        box = "Wm5::WMBox3d"

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
        x.add_method("intersect", "bool", [constrefparam(box, "box")], is_const=True)
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
        x.add_method("distance", "double", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("closestPoint", "Vector3d", [constrefparam("Vector3d", "point")], is_const=True)
        x.add_method("reconstruct", None, [constrefparam("vector_of_Vector3d", "vertices"),
                                           constrefparam("vector_of_vector_of_unsigned", "facetVertices"),
                                           constrefparam("vector_of_Vector3d", "facetNormals")])
        x.add_method("convex", "bool", [param("double", "tol", default_value="1.0e-8")], is_const=True)

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

        return

    #-------------------------------------------------------------------------------
    # Add the free functions for a given dimensionality.
    #-------------------------------------------------------------------------------
    def addDimFunctions(self, space, ndim):

        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim

        space.add_function("computeEigenValues", None, 
                           [constrefparam(symtensorfield, "field"),
                            refparam(vectorfield, "eigenValues"),
                            refparam(tensorfield, "eigenVectors")],
                           template_parameters = ["Dim<%i> " % ndim],
                           custom_name = "computeEigenValues",
                           docstring = "Compute the eigen values & vectors for a field of symmetric tensors.")
        return
