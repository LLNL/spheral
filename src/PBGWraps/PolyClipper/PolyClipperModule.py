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

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addPolygonMethods(self.Polygon)

        self.space.add_function("convertToPolygon", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 constrefparam("Spheral::Polygon", "Spheral_polygon")],
                                docstring = "Convert a Spheral polygon to a PolyClipper polygon.")

        self.space.add_function("convertFromPolygon", None,
                                [refparam("Spheral::Polygon", "Spheral_polygon"),
                                 constrefparam("PolyClipper::Polygon", "polygon")],
                                docstring = "Convert a PolyClipper polygon to a Spheral polygon.")

        self.space.add_function("copyPolygon", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 constrefparam("PolyClipper::Polygon", "polygon0")],
                                docstring = "Copy a PolyClipper polygon to another.")

        self.space.add_function("moments", None,
                                [refparam("double", "zerothMoment"),
                                 refparam("Spheral::Vector2d", "firstMoment"),
                                 constrefparam("PolyClipper::Polygon", "polygon")],
                                docstring = "Compute the zeroth and first moment of a PolyClipper polygon.")

        self.space.add_function("clipPolygon", None,
                                [refparam("PolyClipper::Polygon", "polygon"),
                                 constrefparam("vector_of_Plane2d", "planes")],
                                docstring = "Clip a PolyClipper polygon by a set of planes.")

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
    
        me = "Polygon"

        # Constructors.
        x.add_constructor([])

        return
