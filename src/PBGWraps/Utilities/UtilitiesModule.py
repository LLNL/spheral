from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Utilities:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Utilities/UtilitiesTypes.hh"')
        self.Spheral = mod.add_cpp_namespace("Spheral")
        self.NodeSpace = self.Spheral.add_cpp_namespace("NodeSpace")

        # Expose types.
        self.NewtonRaphsonFunction = addObject(mod, "NewtonRaphsonFunction", allow_subclassing=True)
        self.SimpsonsIntegrationDoubleFunction = addObject(mod, "SimpsonsIntegrationDoubleFunction", allow_subclassing=True)
        self.KeyTraits = addObject(mod, "KeyTraits", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        Spheral = mod.add_cpp_namespace("Spheral")

        self.NewtonRaphsonFunction.add_constructor([])
        self.NewtonRaphsonFunction.add_method("__call__", "pair_double_double",
                                              [param("double", "x")],
                                              is_const=True,
                                              is_pure_virtual=True)

        self.SimpsonsIntegrationDoubleFunction.add_constructor([])
        self.SimpsonsIntegrationDoubleFunction.add_method("__call__", "double",
                                                          [param("double", "x")],
                                                          is_const=True,
                                                          is_pure_virtual=True)

        Spheral.add_function("erff", "double", [param("double", "x")], docstring="You know, the error function.")
        Spheral.add_function("newtonRaphsonFindRoot", "double", [constrefparam("NewtonRaphsonFunction", "function"),
                                                                 param("float", "x1"),
                                                                 param("float", "x2"),
                                                                 param("float", "xaccuracy", default_value="1.0e-15"),
                                                                 param("float", "yaccuracy", default_value="1.0e-15"),
                                                                 param("int", "maxIterations", default_value="100")],
                             docstring="Newton-Raphson root finder.")
        Spheral.add_function("simpsonsIntegrationDouble", "double", [constrefparam("SimpsonsIntegrationDoubleFunction", "function"),
                                                                     param("double", "x0"),
                                                                     param("double", "x1"),
                                                                     param("unsigned int", "numBins")],
                             docstring="Simpsons rule integration for double functions.")

        # Pack/unpack types to strings.
        for value, name in (("int",                          "Int"),
                            ("double",                       "Double"),
                            ("uint32_t",                     "UL"),
                            ("uint64_t",                     "ULL"),
                            ("vector_of_unsigned",           "VectorOfUnsigned"),
                            ("vector_of_int",                "VectorOfInt"),
                            ("vector_of_float",              "VectorOfFloat"),
                            ("vector_of_double",             "VectorOfDouble"),
                            ("vector_of_ULL",                "VectorOfULL"),
                            ("vector_of_vector_of_unsigned", "VectorOfVectorOfUnsigned"),
                            ("Box1d",                        "Box1d"),
                            ("Polygon",                      "Polygon"),
                            ("Polyhedron",                   "Polyhedron"),
                            ):
            exec("""
Spheral.add_function("convertElementToString", "std::string", [param("%(value)s", "x")],
                     template_parameters = ["%(value)s"],
                     custom_name = "packElement%(name)s")
Spheral.add_function("convertStringToElement", "%(value)s", [param("std::string", "x")],
                     template_parameters = ["%(value)s"],
                     custom_name = "unpackElement%(name)s")
""" % {"value" : value, "name" : name})

        # boundingVolumes
        for dim, value, vector, inst in (("1d", "Box1d", "Vector1d", "Dim<1>"),
                                         ("2d", "Polygon", "Vector2d", "Dim<2>"),
                                         ("3d", "Polyhedron", "Vector3d", "Dim<3>")):
            exec("""
Spheral.add_function("boundingBox", None, [constrefparam("vector_of_%(vector)s", "positions"),
                                           refparam("%(vector)s", "xmin"),
                                           refparam("%(vector)s", "xmax")],
                                           template_parameters = ["%(vector)s"],
                                           custom_name = "boundingBox")
Spheral.add_function("globalBoundingBox", None, [constrefparam("Spheral::FieldSpace::VectorFieldList%(dim)s", "positions"),
                                                 refparam("%(vector)s", "xmin"),
                                                 refparam("%(vector)s", "xmax"),
                                                 param("bool", "ghost", default_value="false")],
                                                 template_parameters = ["%(inst)s"],
                                                 custom_name = "globalBoundingBox")
Spheral.add_function("globalBoundingVolumes", None, [constrefparam("DataBase%(dim)s", "dataBase"),
                                                     refparam("%(value)s", "nodeVolume"),
                                                     refparam("%(value)s", "sampleVolume")],
                                                     template_parameters = ["%(inst)s"],
                                                     custom_name = "globalBoundingVolumes")
""" % {"dim" : dim, "value" : value, "vector" : vector, "inst" : inst})

        # Stuff that depends on dimension.
        for dim in ("1d", "2d", "3d"):
            vector = "Vector%s" % dim
            Spheral.add_function("collinear", "bool",
                                 [constrefparam(vector, "a"), constrefparam(vector, "b"), constrefparam(vector, "c"),
                                  param("double", "tol", default_value="1.0e-10")],
                                 template_parameters = [vector],
                                 custom_name = "collinear",
                                 docstring = "Test if the three points are collinear.")
            Spheral.add_function("between", "bool",
                                 [constrefparam(vector, "a"), constrefparam(vector, "b"), constrefparam(vector, "c"), param("double", "tol", default_value="1.0e-10")],
                                 template_parameters = [vector],
                                 custom_name = "between",
                                 docstring = "Test if point c is between a & b.")
            Spheral.add_function("segmentSegmentDistance", "double", [constrefparam(vector, "a0"),
                                                                      constrefparam(vector, "a1"),
                                                                      constrefparam(vector, "b0"),
                                                                      constrefparam(vector, "b1")],
                                 docstring = "Find the distance between line segements (a0,a1) -> (b0,b1)")
            Spheral.add_function("segmentSegmentIntersection", "char", [constrefparam(vector, "a0"),
                                                                        constrefparam(vector, "a1"),
                                                                        constrefparam(vector, "b0"),
                                                                        constrefparam(vector, "b1"),
                                                                        refparam(vector, "result1"),
                                                                        refparam(vector, "result2"),
                                                                        param("double", "tol", default_value="1.0e-8")],
                                 docstring = "Compute the intersection of two line segments (a0,a1) (b0,b1).")
            Spheral.add_function("segmentSegmentIntersection", "bool",
                                 [constrefparam(vector, "a0"),
                                  constrefparam(vector, "a1"),
                                  constrefparam(vector, "b0"),
                                  constrefparam(vector, "b1"),
                                  param("double", "tol", default_value="1.0e-10")],
                                 template_parameters = [vector],
                                 custom_name = "segmentSegmentIntersection"),

        # These methods are only valid in 2d and 3d.
        for dim in ("2d", "3d"):
            vector = "Vector%s" % dim
            Spheral.add_function("pointPlaneDistance", "double", 
                                 [constrefparam(vector, "point"),
                                  constrefparam(vector, "origin"),
                                  constrefparam(vector, "unitNormal")],
                                 template_parameters = [vector],
                                 custom_name = "pointPlaneDistance",
                                 docstring="Compute the distance from (point) to the plane defined by (origin, unitNormal).")
            Spheral.add_function("closestPointOnSegment", vector,
                                 [constrefparam(vector, "p"), constrefparam(vector, "a0"), constrefparam(vector, "a1")],
                                 template_parameters = [vector],
                                 custom_name = "closestPointOnSegment",
                                 docstring = "Find the closest point on a line segment (a0,a1) to point (p).")

        # Closest point in plane to a point.
        Spheral.add_function("closestPointOnPlane", "Vector3d",
                             [constrefparam("Vector3d", "p"), constrefparam("Vector3d", "origin"), constrefparam("Vector3d", "unitNormal")],
                             docstring = "Find the closest point in the plane (origin,normal) to point (p).")

        # Line segment-plane intersections.
        Spheral.add_function("segmentPlaneIntersection", "char", [constrefparam("Vector3d", "s0"),
                                                                  constrefparam("Vector3d", "s1"),
                                                                  constrefparam("Vector3d", "point"),
                                                                  constrefparam("Vector3d", "normal"),
                                                                  refparam("Vector3d", "result"),
                                                                  param("double", "tol", default_value="1.0e-8")],
                             docstring = "Compute the intesection of a line segment (s0,s1) with a plane (point,normal).")
        Spheral.add_function("segmentPlaneIntersection", "char", [constrefparam("Vector3d", "s0"),
                                                                  constrefparam("Vector3d", "s1"),
                                                                  constrefparam("Vector3d", "p0"),
                                                                  constrefparam("Vector3d", "p1"),
                                                                  constrefparam("Vector3d", "p2"),
                                                                  refparam("Vector3d", "result"),
                                                                  param("double", "tol", default_value="1.0e-8")],
                             docstring = "Compute the intesection of a line segment (s0,s1) with a plane (p0,p1,p2).")
        Spheral.add_function("segmentPlanarSectionIntersection", "char", [constrefparam("Vector3d", "s0"),
                                                                          constrefparam("Vector3d", "s1"),
                                                                          constrefparam("vector_of_Vector3d", "pverts"),
                                                                          refparam("Vector3d", "result"),
                                                                          param("double", "tol", default_value="1.0e-8")],
                             docstring = "Compute the intesection of a line segment (s0,s1) with a polygonal planar section (pverts).")
        Spheral.add_function("segmentPlaneIntersection", "bool", [constrefparam("Vector3d", "a0"),
                                                                  constrefparam("Vector3d", "a1"),
                                                                  constrefparam("vector_of_Vector3d", "vertices"),
                                                                  constrefparam("Vector3d", "normal"),
                                                                  param("double", "tol", default_value="1.0e-10")],
                             docstring = "Test if a line segment intersects a planar section.")
        Spheral.add_function("segmentPlaneIntersection", "bool", [constrefparam("Vector3d", "a0"),
                                                                  constrefparam("Vector3d", "a1"),
                                                                  constrefparam("vector_of_Vector3d", "vertices"),
                                                                  constrefparam("vector_of_unsigned", "ipoints"),
                                                                  constrefparam("Vector3d", "normal"),
                                                                  param("double", "tol", default_value="1.0e-10")],
                             docstring = "Test if a line segment intersects a planar section.")

        # Rectlinear mesh stuff.
        Spheral.add_function("writeRectilinearMesh", None, [constrefparam("std::string", "fileName"),
                                                            param("bool", "binaryFile"),
                                                            constrefparam("vector_of_int", "dimensions"),
                                                            constrefparam("vector_of_vector_of_double", "coords"),
                                                            constrefparam("vector_of_string", "scalarNames"),
                                                            constrefparam("vector_of_string", "vectorNames"),
                                                            constrefparam("vector_of_string", "tensorNames"),
                                                            constrefparam("vector_of_string", "symTensorNames"),
                                                            constrefparam("vector_of_vector_of_double", "sampledScalars"),
                                                            constrefparam("vector_of_vector_of_Vector2d", "sampledVectors"),
                                                            constrefparam("vector_of_vector_of_Tensor2d", "sampledTensors"),
                                                            constrefparam("vector_of_vector_of_SymTensor2d", "sampledSymTensors")],
                             template_parameters = ["Dim<2>"],
                             custom_name = "writeRectilinearMesh2d",
                             docstring="Write out an XY VISIT rectilinear mesh file.")
        Spheral.add_function("writeRectilinearMesh", None, [constrefparam("std::string", "fileName"),
                                                            param("bool", "binaryFile"),
                                                            constrefparam("vector_of_int", "dimensions"),
                                                            constrefparam("vector_of_vector_of_double", "coords"),
                                                            constrefparam("vector_of_string", "scalarNames"),
                                                            constrefparam("vector_of_string", "vectorNames"),
                                                            constrefparam("vector_of_string", "tensorNames"),
                                                            constrefparam("vector_of_string", "symTensorNames"),
                                                            constrefparam("vector_of_vector_of_double", "sampledScalars"),
                                                            constrefparam("vector_of_vector_of_Vector3d", "sampledVectors"),
                                                            constrefparam("vector_of_vector_of_Tensor3d", "sampledTensors"),
                                                            constrefparam("vector_of_vector_of_SymTensor3d", "sampledSymTensors")],
                             template_parameters = ["Dim<3>"],
                             custom_name = "writeRectilinearMesh3d",
                             docstring="Write out an XYZ VISIT rectilinear mesh file.")
                                                            

        # Polygon/Polyhedron containment.
        for volume, name, vector in (("Polygon", "polygon", "Vector2d"),
                                     ("Polyhedron", "polyhedron", "Vector3d")):
            exec("""
Spheral.add_function("pointOn%(volume)s", "bool", [constrefparam("%(vector)s", "p"),
                                                   constrefparam("vector_of_%(vector)s", "vertices"),
                                                   param("double", "tol", default_value="1.0e-10")],
                     docstring = "Test if the given point is on the boundary of a %(name)s specified by it's vertices.")
Spheral.add_function("pointIn%(volume)s", "bool", [constrefparam("%(vector)s", "p"),
                                                   constrefparam("vector_of_%(vector)s", "vertices"),
                                                   param("bool", "countBoundary", default_value="false"),
                                                   param("double", "tol", default_value="1.0e-10")],
                     docstring = "Test if the given point is contained withing a %(name)s specified by it's vertices.")
Spheral.add_function("pointIn%(volume)s", "bool", [constrefparam("%(vector)s", "p"),
                                                   constrefparam("%(volume)s", "%(name)s"),
                                                   param("bool", "countBoundary", default_value="false"),
                                                   param("double", "tol", default_value="1.0e-10")],
                     docstring = "Test if the given point is contained withing a %(name)s.")
Spheral.add_function("segmentIntersectEdges", "bool", [constrefparam("%(vector)s", "a0"),
                                                       constrefparam("%(vector)s", "a1"),
                                                       constrefparam("%(volume)s", "poly"),
                                                       param("double", "tol", default_value="1.0e-8")],
                     docstring = "Test if the give line segment intersects any edges/vertices of the given faceted volume.")
""" % {"volume" : volume, "name" : name, "vector" : vector})

        Spheral.add_function("pointInPolygon", "bool", 
                             [constrefparam("Vector3d", "p"), 
                              constrefparam("vector_of_Vector3d", "vertices"),
                              constrefparam("vector_of_unsigned", "ipoints"),
                              constrefparam("Vector3d", "normal")],
                             docstring = "Test if the given 3-D point p is in a polygon.")

        Spheral.add_function("refinePolyhedron", "Polyhedron", 
                             [constrefparam("Polyhedron", "poly0"), param("int", "numLevels")],
                             docstring = "Return a new Polyhedron based on refining an existing one a given number of levels.")

        # Boost.math functions.
        Spheral.add_function("legendre_p", "double", 
                             [param("int", "l"), param("int", "m"), param("double", "x")],
                             docstring = "Compute the associated Legendre polynomial.")

        # Add the KeyTraits attributes.
        self.KeyTraits.add_static_attribute("numbits", "int",  is_const=True)
        self.KeyTraits.add_static_attribute("numbits1d", "int",  is_const=True)
        self.KeyTraits.add_static_attribute("zero", "uint64_t",  is_const=True)
        self.KeyTraits.add_static_attribute("one", "uint64_t",  is_const=True)
        self.KeyTraits.add_static_attribute("two", "uint64_t",  is_const=True)
        self.KeyTraits.add_static_attribute("maxKey1d", "uint64_t",  is_const=True)
        self.KeyTraits.add_static_attribute("maxKey", "uint64_t",  is_const=True)

        # Dimension dependent bindings.
        self.generateDimBindings(1)
        self.generateDimBindings(2)
        self.generateDimBindings(3)

        return

    #---------------------------------------------------------------------------
    # Generate dimension dependent bindings.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):

        dim = "Dim<%i> " % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        nodelist = "NodeList%id" % ndim
        database = "DataBase%id" % ndim
        intfield = "IntField%id" % ndim
        intfieldlist = "IntFieldList%id" % ndim
        ullfieldlist = "ULLFieldList%id" % ndim
        scalarfieldlist = "ScalarFieldList%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "TableKernel%id" % ndim
        smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim

        Spheral = self.Spheral
        NodeSpace = self.NodeSpace

        # Expose methods.
        NodeSpace.add_function("numGlobalNodes", "int", [constrefparam(nodelist, "nodes")],
                               custom_name = "numGlobalNodes%id" % ndim,
                               template_parameters = [dim],
                               docstring="Global number of nodes in the NodeList.")

        NodeSpace.add_function("numGlobalNodesAll", "int", [constrefparam(database, "dataBase")],
                               custom_name = "numGlobalNodesAll%id" % ndim,
                               template_parameters = [dim],
                               docstring="Global number of nodes in the DataBase.")

        NodeSpace.add_function("globalNodeIDs", intfield, [constrefparam(nodelist, "nodes")],
                               custom_name = "globalNodeIDs%id" % ndim,
                               template_parameters = [dim],
                               docstring="Determine unique global node IDs for the nodes in a NodeList.")

        NodeSpace.add_function("globalNodeIDsAll", intfieldlist, [constrefparam(database, "dataBase")],
                               custom_name = "globalNodeIDsAll%id" % ndim,
                               template_parameters = [dim],
                               docstring="Determine unique global node IDs for the nodes in a DataBase.")

        Spheral.add_function("rotationMatrix%id" % ndim, tensor, [constrefparam(vector, "runit")],
                             docstring="Rotational transformation to align with the given unit vector.")

        Spheral.add_function("iterateIdealH", None, [constrefparam(database, "dataBase"),
                                                     constrefparam(vector_of_boundary, "boundaries"),
                                                     constrefparam(tablekernel, "W"),
                                                     constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                                                     param("int", "maxIterations", default_value="100"),
                                                     param("double", "tolerance", default_value="1.0e-10"),
                                                     param("double", "nPerhForIteration", default_value="0.0"),
                                                     param("int", "sphericalStart", default_value="0"),
                                                     param("int", "fixDeterminant", default_value="0")],
                             template_parameters = [dim],
                             custom_name = "iterateIdealH%id" % ndim,
                             docstring="Iterate the Hfield for NodeLists in the DataBase using the ideal H algorithm.")

        Spheral.add_function("mortonOrderIndicies%id" % ndim, ullfieldlist,
                             [constrefparam(database, "dataBase")],
                             docstring = "Compute indicies for nodes obeying Morton ordering.")
        
        Spheral.add_function("peanoHilbertOrderIndicies%id" % ndim, ullfieldlist,
                             [constrefparam(database, "dataBase")],
                             docstring = "Compute indicies for nodes obeying the Peano-Hilbert ordering.")
        

        Spheral.add_function("integrateThroughMeshAlongSegment%id" % ndim,
                             "double",
                             [constrefparam("vector_of_vector_of_double", "values"),
                              constrefparam(vector, "xmin"),
                              constrefparam(vector, "xmax"),
                              constrefparam("vector_of_unsigned", "ncells"),
                              constrefparam(vector, "s0"),
                              constrefparam(vector, "s1")],
                             docstring = "Integrate through a lattice sampled field along a line segment.")

        Spheral.add_function("numberDensity",
                             scalarfieldlist,
                             [constrefparam(database, "dataBase"), constrefparam(tablekernel, "W")],
                             template_parameters = [dim],
                             custom_name = "numberDensity%id" % ndim,
                             docstring = "Compute the ASPH sum number density for each node in a DataBase.")

        Spheral.add_function("testBoxIntersection%id" % ndim,
                             "bool",
                             [constrefparam(vector, "xmin1"), constrefparam(vector, "xmax1"),
                              constrefparam(vector, "xmin2"), constrefparam(vector, "xmax2")],
                             docstring = "Test if the two boxes intersect.")
                              
        Spheral.add_function("testPointInBox%id" % ndim,
                             "bool",
                             [constrefparam(vector, "point"),
                              constrefparam(vector, "xmin"), constrefparam(vector, "xmax")],
                             docstring = "Test if the point is in the box.")
                              
        Spheral.add_function("planarReflectingOperator",
                             tensor,
                             [constrefparam(plane, "plane")],
                             template_parameters = [plane],
                             custom_name = "planarReflectingOperator%id" % ndim,
                             docstring = "Generate the planar reflection transformation for th given plane.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

