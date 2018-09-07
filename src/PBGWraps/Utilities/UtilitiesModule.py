from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Utilities:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/UtilitiesTypes.hh"' % srcdir)
        self.Spheral = mod.add_cpp_namespace("Spheral")
        self.PythonBoundFunctors = self.Spheral.add_cpp_namespace("PythonBoundFunctors")

        # Expose types.
        self.KeyTraits = addObject(mod, "KeyTraits", allow_subclassing=True)
        self.Timer = addObject(mod, "Timer")
        self.ScalarScalarFunctor = addObject(self.PythonBoundFunctors, "SpheralFunctor", template_parameters=["double", "double"], custom_name="ScalarScalarFunctor", allow_subclassing=True)
        self.ScalarPairScalarFunctor = addObject(self.PythonBoundFunctors, "SpheralFunctor", template_parameters=["double", "std::pair<double, double>"], custom_name="ScalarPairScalarFunctor", allow_subclassing=True)
        for ndim in self.dims:
            exec("""
self.VectorScalarFunctor%(ndim)id = addObject(self.PythonBoundFunctors, "SpheralFunctor", template_parameters=["Vector%(ndim)id", "double"], custom_name="VectorScalarFunctor%(ndim)id", allow_subclassing=True)
self.VectorVectorFunctor%(ndim)id = addObject(self.PythonBoundFunctors, "SpheralFunctor", template_parameters=["Vector%(ndim)id", "Vector%(ndim)id"], custom_name="VectorVectorFunctor%(ndim)id", allow_subclassing=True)
self.VectorPairScalarFunctor%(ndim)id = addObject(self.PythonBoundFunctors, "SpheralFunctor", template_parameters=["Vector%(ndim)id", "std::pair<double, double>"], custom_name="VectorPairScalarFunctor%(ndim)id", allow_subclassing=True)
""" % {"ndim" : ndim})

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        Spheral = mod.add_cpp_namespace("Spheral")

        self.addTimerBindings(self.Timer)

        # Add the functors.
        self.addFunctorBindings(self.ScalarScalarFunctor,       "double",   "double")
        self.addFunctorBindings(self.ScalarPairScalarFunctor,   "double",   "pair_double_double")

        Spheral.add_function("erff", "double", [param("double", "x")], docstring="You know, the error function.")
        Spheral.add_function("newtonRaphsonFindRoot", "double", [constrefparam("Spheral::PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double> >", "function"),
                                                                 param("float", "x1"),
                                                                 param("float", "x2"),
                                                                 param("float", "xaccuracy", default_value="1.0e-15"),
                                                                 param("float", "yaccuracy", default_value="1.0e-15"),
                                                                 param("int", "maxIterations", default_value="100")],
                             docstring="Newton-Raphson root finder.")
        Spheral.add_function("simpsonsIntegrationDouble", "double", [constrefparam("Spheral::PythonBoundFunctors::SpheralFunctor<double, double>", "function"),
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
        for ndim in self.dims:
            exec("""
Spheral.add_function("boundingBox", None, [constrefparam("vector_of_%(vector)s", "positions"),
                                           refparam("%(vector)s", "xmin"),
                                           refparam("%(vector)s", "xmax")],
                                           template_parameters = ["%(vector)s"],
                                           custom_name = "boundingBox")
Spheral.add_function("globalBoundingBox", None, [constrefparam("Spheral::VectorFieldList%(dim)s", "positions"),
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
""" % {"dim"   : "%id" % ndim, 
       "value" : ("Box1d", "Polygon", "Polyhedron")[ndim - 1],
       "vector": "Vector%id" % ndim,
       "inst"  : "Dim<%i>" % ndim})

        # Stuff that depends on dimension.
        for ndim in self.dims:
            dim = "%id" % ndim
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
        for ndim in (2, 3):
            if ndim in self.dims:
                dim = "Dim<%i>" % ndim
                vector = "Vector%id" % ndim
                vector_of_boundary = "vector_of_Boundary%id" % ndim
                vector_of_scalarfieldptr = "vector_of_ScalarFieldPtr%id" % ndim
                vector_of_vectorfieldptr = "vector_of_VectorFieldPtr%id" % ndim
                vector_of_tensorfieldptr = "vector_of_TensorFieldPtr%id" % ndim
                vector_of_symTensorfieldptr = "vector_of_SymTensorFieldPtr%id" % ndim
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
                Spheral.add_function("overlayRemapFields", None,
                                     [constrefparam(vector_of_boundary, "boundaries"),
                                      constrefparam(vector_of_scalarfieldptr, "scalarDonorFields"),
                                      constrefparam(vector_of_vectorfieldptr, "vectorDonorFields"),
                                      constrefparam(vector_of_tensorfieldptr, "tensorDonorFields"),
                                      constrefparam(vector_of_symTensorfieldptr, "symTensorDonorFields"),
                                      refparam(vector_of_scalarfieldptr, "scalarAcceptorFields"),
                                      refparam(vector_of_vectorfieldptr, "vectorAcceptorFields"),
                                      refparam(vector_of_tensorfieldptr, "tensorAcceptorFields"),
                                      refparam(vector_of_symTensorfieldptr, "symTensorAcceptorFields")],
                                     template_parameters = [dim],
                                     custom_name = "overlayRemapFields",
                                     docstring = "Do a simple donor overlay using geometric intersection.")

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

        # # Rectlinear mesh stuff.
        # Spheral.add_function("writeRectilinearMesh", None, [constrefparam("std::string", "fileName"),
        #                                                     param("bool", "binaryFile"),
        #                                                     constrefparam("vector_of_int", "dimensions"),
        #                                                     constrefparam("vector_of_vector_of_double", "coords"),
        #                                                     constrefparam("vector_of_string", "scalarNames"),
        #                                                     constrefparam("vector_of_string", "vectorNames"),
        #                                                     constrefparam("vector_of_string", "tensorNames"),
        #                                                     constrefparam("vector_of_string", "symTensorNames"),
        #                                                     constrefparam("vector_of_vector_of_double", "sampledScalars"),
        #                                                     constrefparam("vector_of_vector_of_Vector2d", "sampledVectors"),
        #                                                     constrefparam("vector_of_vector_of_Tensor2d", "sampledTensors"),
        #                                                     constrefparam("vector_of_vector_of_SymTensor2d", "sampledSymTensors")],
        #                      template_parameters = ["Dim<2>"],
        #                      custom_name = "writeRectilinearMesh2d",
        #                      docstring="Write out an XY VISIT rectilinear mesh file.")
        # Spheral.add_function("writeRectilinearMesh", None, [constrefparam("std::string", "fileName"),
        #                                                     param("bool", "binaryFile"),
        #                                                     constrefparam("vector_of_int", "dimensions"),
        #                                                     constrefparam("vector_of_vector_of_double", "coords"),
        #                                                     constrefparam("vector_of_string", "scalarNames"),
        #                                                     constrefparam("vector_of_string", "vectorNames"),
        #                                                     constrefparam("vector_of_string", "tensorNames"),
        #                                                     constrefparam("vector_of_string", "symTensorNames"),
        #                                                     constrefparam("vector_of_vector_of_double", "sampledScalars"),
        #                                                     constrefparam("vector_of_vector_of_Vector3d", "sampledVectors"),
        #                                                     constrefparam("vector_of_vector_of_Tensor3d", "sampledTensors"),
        #                                                     constrefparam("vector_of_vector_of_SymTensor3d", "sampledSymTensors")],
        #                      template_parameters = ["Dim<3>"],
        #                      custom_name = "writeRectilinearMesh3d",
        #                      docstring="Write out an XYZ VISIT rectilinear mesh file.")
                                                            

        # Polygon/Polyhedron containment.
        for volume, name, vector in (("Polygon", "polygon", "Vector2d"),
                                     ("Polyhedron", "polyhedron", "Vector3d")):
            exec("""
Spheral.add_function("pointOn%(volume)s", "bool", [constrefparam("%(vector)s", "p"),
                                                   constrefparam("%(volume)s", "%(name)s"),
                                                   param("double", "tol", default_value="1.0e-10")],
                     docstring = "Test if the given point is on the boundary of a %(name)s specified by it's vertices.")
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

        # R2D/R3D utilities
        Spheral.add_function("clipFacetedVolume", "Polygon", 
                             [constrefparam("Polygon", "poly"),
                              constrefparam("vector_of_Plane2d", "planes")],
                             docstring = "Clip a polygon with a set of planes.")
        Spheral.add_function("clipFacetedVolume", "Polyhedron", 
                             [constrefparam("Polyhedron", "poly"),
                              constrefparam("vector_of_Plane3d", "planes")],
                             docstring = "Clip a polyhedron with a set of planes.")
        Spheral.add_function("clippedVolume", "double", 
                             [constrefparam("Polygon", "poly"),
                              constrefparam("vector_of_Plane2d", "planes")],
                             docstring = "Return the area of a polygon clipped with a set of planes.")
        Spheral.add_function("clippedVolume", "double", 
                             [constrefparam("Polyhedron", "poly"),
                              constrefparam("vector_of_Plane3d", "planes")],
                             docstring = "Return the volume of a polyhedron clipped with a set of planes.")

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
        for dim in self.dims:
            self.generateDimBindings(dim)
        for dim in (1, 2, 3):
            self.generateAllDimBindings(dim)

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
        vectorfieldlist = "VectorFieldList%id" % ndim
        tensorfieldlist = "TensorFieldList%id" % ndim
        symtensorfieldlist = "SymTensorFieldList%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "TableKernel%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        Spheral = self.Spheral

        # Dimension dependent functor bindings.
        self.addFunctorBindings(eval("self.VectorScalarFunctor%id" % ndim), "Vector%id" % ndim, "double")
        self.addFunctorBindings(eval("self.VectorVectorFunctor%id" % ndim), "Vector%id" % ndim, "Vector%id" % ndim)
        self.addFunctorBindings(eval("self.VectorPairScalarFunctor%id" % ndim), "Vector%id" % ndim, "pair_double_double")

        # Expose methods.
        Spheral.add_function("numGlobalNodes", "int", [constrefparam(nodelist, "nodes")],
                               custom_name = "numGlobalNodes%id" % ndim,
                               template_parameters = [dim],
                               docstring="Global number of nodes in the NodeList.")

        Spheral.add_function("numGlobalNodesAll", "int", [constrefparam(database, "dataBase")],
                               custom_name = "numGlobalNodesAll%id" % ndim,
                               template_parameters = [dim],
                               docstring="Global number of nodes in the DataBase.")

        Spheral.add_function("globalNodeIDs", intfield, [constrefparam(nodelist, "nodes")],
                               custom_name = "globalNodeIDs%id" % ndim,
                               template_parameters = [dim],
                               docstring="Determine unique global node IDs for the nodes in a NodeList.")

        Spheral.add_function("globalNodeIDsAll", intfieldlist, [constrefparam(database, "dataBase")],
                               custom_name = "globalNodeIDsAll%id" % ndim,
                               template_parameters = [dim],
                               docstring="Determine unique global node IDs for the nodes in a DataBase.")

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

        Spheral.add_function("nodeOrdering%id" % ndim, intfieldlist,
                             [constrefparam(ullfieldlist, "criteria")],
                             docstring = "Return a node ordering index list according to the given comparison criteria.")

        Spheral.add_function("mortonOrderIndices%id" % ndim, ullfieldlist,
                             [constrefparam(vectorfieldlist, "positions")],
                             docstring = "Compute indices for nodes obeying Morton ordering.")
        
        Spheral.add_function("mortonOrderIndices%id" % ndim, ullfieldlist,
                             [constrefparam(database, "dataBase")],
                             docstring = "Compute indices for nodes obeying Morton ordering.")
        
        Spheral.add_function("peanoHilbertOrderIndices%id" % ndim, ullfieldlist,
                             [constrefparam(vectorfieldlist, "positions")],
                             docstring = "Compute indices for nodes obeying the Peano-Hilbert ordering.")

        Spheral.add_function("peanoHilbertOrderIndices%id" % ndim, ullfieldlist,
                             [constrefparam(database, "dataBase")],
                             docstring = "Compute indices for nodes obeying the Peano-Hilbert ordering.")

        Spheral.add_function("numberDensity",
                             scalarfieldlist,
                             [constrefparam(database, "dataBase"), constrefparam(tablekernel, "W")],
                             template_parameters = [dim],
                             custom_name = "numberDensity%id" % ndim,
                             docstring = "Compute the ASPH sum number density for each node in a DataBase.")

        Spheral.add_function("integrateThroughMeshAlongSegment%id" % ndim,
                             "double",
                             [constrefparam("vector_of_vector_of_double", "values"),
                              constrefparam(vector, "xmin"),
                              constrefparam(vector, "xmax"),
                              constrefparam("vector_of_unsigned", "ncells"),
                              constrefparam(vector, "s0"),
                              constrefparam(vector, "s1")],
                             docstring = "Integrate through a lattice sampled field along a line segment.")

        for fl in (scalarfieldlist, vectorfieldlist, tensorfieldlist, symtensorfieldlist):
            Spheral.add_function("computeShepardsInterpolation",
                                 fl,
                                 [constrefparam(fl, "fieldList"),
                                  constrefparam(connectivitymap, "connectivityMap"),
                                  constrefparam(tablekernel, "W"),
                                  constrefparam(vectorfieldlist, "position"),
                                  constrefparam(symtensorfieldlist, "H"),
                                  constrefparam(scalarfieldlist, "weight")],
                                 docstring = "Interpolate a FieldList using a Shepards function approach.")

        return

    #---------------------------------------------------------------------------
    # Generate dimension dependent bindings.
    #---------------------------------------------------------------------------
    def generateAllDimBindings(self, ndim):

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
        vectorfieldlist = "VectorFieldList%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "TableKernel%id" % ndim
        smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim

        Spheral = self.Spheral

        # Expose methods.
        Spheral.add_function("rotationMatrix%id" % ndim, tensor, [constrefparam(vector, "runit")],
                             docstring="Rotational transformation to align with the given unit vector.")

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
    # Timer bindings.
    #---------------------------------------------------------------------------
    def addTimerBindings(self, x):

        # Constructors
        x.add_constructor([])
        x.add_constructor([param("std::string", "name")])
        x.add_constructor([param("std::string", "name"),
                           refparam("Timer", "base")])

        # Methods
        x.add_method("setup", None, [])
        x.add_method("start", None, [])
        x.add_method("stop", None, [])
        x.add_method("clear", None, [])
        x.add_method("getTimeStampWC", "double", [])
        x.add_method("wc_time", "double", [])
        x.add_method("Name", "std::string", [])
        x.add_method("Count", "long int", [])

        # Static methods
        x.add_method("TimerSummary", None, [], is_static=True)

        return

    #---------------------------------------------------------------------------
    # functor bindings.
    #---------------------------------------------------------------------------
    def addFunctorBindings(self, F, argval, retval):
        F.add_constructor([])
        F.add_method("__call__", retval,
                     [param(argval, "x")],
                     is_const=True,
                     is_pure_virtual=True)
        return
