from pybindgen import *

from PBGutils import *
from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class PolyClipper:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"%s/PolyClipperTypes.hh"' % srcdir)

        # Namespace.
        self.space = mod.add_cpp_namespace("PolyClipper")

        # Expose types.
        self.Plane2d = addObject(self.space, "PolyClipperPlane2d")
        self.Plane3d = addObject(self.space, "PolyClipperPlane3d")
        self.Polygon = addObject(self.space, "Polygon")
        self.Polyhedron = addObject(self.space, "Polyhedron")
        self.vector_of_Plane2d = addObject(mod, "vector_of_PolyClipperPlane2d")
        self.vector_of_Plane3d = addObject(mod, "vector_of_PolyClipperPlane3d")

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addPlaneMethods(self.Plane2d, 2)
        self.addPlaneMethods(self.Plane3d, 3)
        self.addPolygonMethods(self.Polygon)
        self.addPolyhedronMethods(self.Polyhedron)

        # vector_of_Plane methods.
        generateStdVectorBindings(self.vector_of_Plane2d, "PolyClipper::PolyClipperPlane2d", "vector_of_PolyClipperPlane2d", indexAsPointer=True)
        generateStdVectorBindings(self.vector_of_Plane3d, "PolyClipper::PolyClipperPlane3d", "vector_of_PolyClipperPlane3d", indexAsPointer=True)

        # Polygon functions
        self.space.add_function("initializePolygon", None,
                                [refparam("PolyClipper::Polygon", "poly"),
                                 constrefparam("vector_of_Vector2d", "position"),
                                 constrefparam("vector_of_vector_of_int", "neighbors")],
                                docstring = "Initialize a polygon from the positions and neighbors of each vertex.")

        self.space.add_function("polygon2string", "std::string",
                                [constrefparam("PolyClipper::Polygon", "polygon")],
                                docstring = "Print a polygon as a formatted string for human reading.")

        self.space.add_function("convertToPolygon", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 constrefparam("Spheral::Polygon", "Spheral_polygon")],
                                docstring = "Convert a Spheral polygon to a PolyClipper polygon.")

        self.space.add_function("convertFromPolygon", None,
                                [refparam("Spheral::Polygon", "Spheral_polygon"),
                                 constrefparam("PolyClipper::Polygon", "polygon")],
                                docstring = "Convert a PolyClipper polygon to a Spheral polygon.")

        self.space.add_function("moments", None,
                                [Parameter.new("double&", "zerothMoment", direction=Parameter.DIRECTION_OUT),
                                 Parameter.new("Spheral::Vector2d&", "zerothMoment", direction=Parameter.DIRECTION_OUT),
                                 constrefparam("PolyClipper::Polygon", "polygon")],
                                docstring = "Compute the zeroth and first moment of a PolyClipper polygon.")

        self.space.add_function("clipPolygon", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 constrefparam("vector_of_PolyClipperPlane2d", "planes")],
                                docstring = "Clip a PolyClipper polygon by a set of planes.")

        self.space.add_function("collapseDegenerates", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 param("double", "tol")],
                                docstring = "Remove degnerate edges/vertices from a PolyClipper polygon.")

        # Polyhedron functions
        self.space.add_function("polyhedron2string", "std::string",
                                [constrefparam("PolyClipper::Polyhedron", "polyhedron")],
                                docstring = "Print a polyhedron as a formatted string for human reading.")

        self.space.add_function("convertToPolyhedron", None,
                                [refparam("PolyClipper::Polyhedron", "polyhedron"),
                                 constrefparam("Spheral::Polyhedron", "Spheral_polyhedron")],
                                docstring = "Convert a Spheral polyhedron to a PolyClipper polyhedron.")

        self.space.add_function("convertFromPolyhedron", None,
                                [refparam("Spheral::Polyhedron", "Spheral_polyhedron"),
                                 constrefparam("PolyClipper::Polyhedron", "polyhedron")],
                                docstring = "Convert a PolyClipper polyhedron to a Spheral polyhedron.")

        self.space.add_function("moments", None,
                                [Parameter.new("double&", "zerothMoment", direction=Parameter.DIRECTION_OUT),
                                 Parameter.new("Spheral::Vector3d&", "zerothMoment", direction=Parameter.DIRECTION_OUT),
                                 constrefparam("PolyClipper::Polyhedron", "polyhedron")],
                                docstring = "Compute the zeroth and first moment of a PolyClipper polyhedron.")

        self.space.add_function("clipPolyhedron", None,
                                [refparam("PolyClipper::Polyhedron", "polyhedron"),
                                 constrefparam("vector_of_PolyClipperPlane3d", "planes")],
                                docstring = "Clip a PolyClipper polyhedron by a set of planes.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["PolyClipper"]

    #-------------------------------------------------------------------------------
    # Plane
    #-------------------------------------------------------------------------------
    def addPlaneMethods(self, x, ndim):
    
        me = "PolyClipperPlane%id" % ndim
        vector = "Spheral::Vector%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("double", "dist"), constrefparam(vector, "normal")])
        x.add_constructor([param(vector, "point"), constrefparam(vector, "normal")])
        x.add_constructor([constrefparam(me, "plane")])

        # Attributes
        x.add_instance_attribute("dist", "double", False)
        x.add_instance_attribute("normal", vector, False)

        return

    #-------------------------------------------------------------------------------
    # Polygon
    #-------------------------------------------------------------------------------
    def addPolygonMethods(self, x):
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam("PolyClipper::Polygon", "poly")])

        # Methods.
        x.add_method("size", "unsigned int", [], is_const=True)

        return

    #-------------------------------------------------------------------------------
    # Polyhedron
    #-------------------------------------------------------------------------------
    def addPolyhedronMethods(self, x):
    
        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam("PolyClipper::Polyhedron", "poly")])

        # Methods.
        x.add_method("size", "unsigned int", [], is_const=True)

        return
