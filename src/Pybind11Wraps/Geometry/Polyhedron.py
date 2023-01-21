#-------------------------------------------------------------------------------
# Polyhedron
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11cppname("GeomPolyhedron")
class Polyhedron:

    PYB11typedefs = """
    using Polyhedron = GeomPolyhedron;
    using Vector = GeomPolyhedron::Vector;
    using Tensor = GeomPolyhedron::Tensor;
    using Facet = GeomPolyhedron::Facet;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit0(self,
                rhs = "const Polyhedron&"):
        "Copy constructor"

    def pyinit1(self,
                points = "const std::vector<Vector>&"):
        """Note this constructor constructs the convex hull of the given points,
meaning that the full set of points passed in may not appear in the vertices."""

    def pyinit2(self,
                points = "const std::vector<Vector>&",
                facetIndices = "const std::vector<std::vector<unsigned> >&"):
        "Construct with explicit vertices and facets"

    @PYB11implementation("[](py::list& points) { std::vector<Vector> vpoints; for (auto p: points) vpoints.push_back(p.cast<Vector>()); return new GeomPolyhedron(vpoints); }")
    def pyinit3(self,
                points = "py::list"):
        "Construct as the convex hull of a python list of points"

    #...........................................................................
    # Methods
    @PYB11const
    def contains(self,
                 point = "const Vector&",
                 countBoundary = ("const bool", "true"),
                 tol = ("const double", "1.0e-8"),
                 useAxom = ("const bool", "false")):
        "Test if the given point is internal to the polyhedron."
        return "bool"

    @PYB11const
    def convexContains(self,
                       point = "const Vector&",
                       countBoundary = ("const bool", "true"),
                       tol = ("const double", "1.0e-8")):
        "Test if the given point is internal to the polyhedron (assumes convexity)."
        return "bool"

    @PYB11const
    def intersect(self,
                  rhs = "const Polyhedron&"):
        "Test if we intersect another polyhedron."
        return "bool"

    @PYB11const
    def convexIntersect(self,
                        rhs = "const Polyhedron&"):
        "Test if we intersect another polyhedron (assumes convexity)"
        return "bool"

    @PYB11const
    @PYB11pycppname("intersect")
    def intersect1(self,
                   rhs = "const std::pair<Vector, Vector>&"):
        "Test if we intersect a box represented by a min/max pair of coordinates."
        return "bool"

    @PYB11const
    @PYB11pycppname("intersect")
    def intersect2(self,
                   s0 = "const Vector&",
                   s1 = "const Vector&"):
        "Test if we intersect a line segment (interior counts as intersection)."
        return "bool"

    @PYB11const
    @PYB11implementation("[](const Polyhedron& self, const Vector& s0, const Vector& s1) { std::vector<unsigned> facetIDs; std::vector<Vector> intersections; self.intersections(s0, s1, facetIDs, intersections); return py::make_tuple(facetIDs, intersections); }")
    def intersections(self,
                      s0 = "const Vector&",
                      s1 = "const Vector&"):
        "Return the intersections of this polyhedron with a line segment denoted by it's end points."
        return "py::tuple"

    @PYB11const
    def edges(self):
        "Get the edges as integer (node) pairs."
        return "std::vector<std::pair<unsigned, unsigned> >"

    def reconstruct(self,
                    vertices = "const std::vector<Vector>&",
                    facetVertices = "const std::vector<std::vector<unsigned> >&"):
        """Reconstruct the internal data given a set of vertices, vertex
indices that define the facets, and outward normals at the facets."""
        return "void"

    @PYB11const
    def closestFacet(self, p = "const Vector&"):
        "Find the facet closest to the given point."
        return "unsigned"

    @PYB11const
    def distance(self, p="const Vector&",
                 useAxom = ("const bool", "false")):
        "Compute the minimum distance to a point.  Note: when useAxom=True this is a signed distance, and otherwise it's positive."
        return "double"

    @PYB11const
    def closestPoint(self, p="const Vector&"):
        "Find the point in the polyhedron closest to the given point."
        return "Vector"

    @PYB11const
    def convex(self, tol = "double"):
        "Test if the polyhedron is convex"
        return "bool"

    def setBoundingBox(self):
        "Set the internal bounding box"
        return "void"

    @PYB11const
    def facetArea(self, facetID="const unsigned"):
        return "double"

    @PYB11const
    def facetAreaNormal(self, facetID="const unsigned"):
        return "Vector"

    @PYB11const
    def facetSubVolume(self, facetID="const unsigned"):
        "Decompose the polyhedron into tetrahedra for each facet"
        return "Polyhedron"

    def transform(self, t="const Tensor&"):
        "Apply a general transformation tensor"
        return "Polyhedron&"

    #...........................................................................
    # Operators
    def __iadd__(self, rhs="Vector()"):
        return

    def __isub__(self, rhs="Vector()"):
        return

    def __add__(self, rhs="Vector()"):
        return

    def __sub__(self, rhs="Vector()"):
        return

    def __imul__(self, rhs="double()"):
        return
    
    def __idiv__(self, rhs="double()"):
        return
    
    def __mul__(self, rhs="double()"):
        return
    
    def __div__(self, rhs="double()"):
        return
    
    def __eq__(self):
        return

    def __ne__(self):
        return

    #...........................................................................
    # Properties
    centroid = PYB11property("Vector")
    vertices = PYB11property("const std::vector<Vector>&", returnpolicy="reference_internal")
    facets = PYB11property("const std::vector<Facet>&", returnpolicy="reference_internal")
    facetVertices = PYB11property("std::vector<std::vector<unsigned> >",
                                  doc="Spit out a vector<vector<unsigned> > that encodes the facets.")
    facetNormals = PYB11property("std::vector<Vector>",
                                 doc="Compute effective surface normals at the facets.")
    facetVertices = PYB11property("std::vector<std::vector<unsigned> >",
                                  doc="Spit out a vector<vector<unsigned> > that encodes the facets.")
    vertexUnitNorms = PYB11property("const std::vector<Vector>&", returnpolicy="reference_internal")
    vertexFacetConnectivity = PYB11property("const std::vector<std::vector<unsigned> >&", returnpolicy="reference_internal")
    facetFacetConnectivity = PYB11property("const std::vector<std::vector<unsigned> >&", returnpolicy="reference_internal")
    xmin = PYB11property("const Vector&", returnpolicy="reference_internal")
    xmax = PYB11property("const Vector&", returnpolicy="reference_internal")
    volume = PYB11property("double")
