from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class PolyClipper:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"Geometry/polyclipper.hh"')

        # Namespace.
        self.space = mod.add_cpp_namespace("PolyClipper")

        # Expose types.
        self.Polygon = addObject(self.space, "Polygon")
        self.Polyhedron = addObject(self.space, "Polyhedron")

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addPolygonMethods(self.Polygon)
        self.addPolyhedronMethods(self.Polyhedron)

        # Polygon functions
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
                                 constrefparam("vector_of_Plane2d", "planes")],
                                docstring = "Clip a PolyClipper polygon by a set of planes.")

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
                                 constrefparam("vector_of_Plane3d", "planes")],
                                docstring = "Clip a PolyClipper polyhedron by a set of planes.")

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["PolyClipper"]

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
