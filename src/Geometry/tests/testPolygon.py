#ATS:test(SELF, label="Polygon unit tests")
# Unit tests for the Polygon class

import unittest
from math import *
from SpheralTestUtilities import fuzzyEqual

from Spheral2d import *

# Create a global random number generator.
import random
rangen = random.Random()

plots = []

#===============================================================================
# Generate random points in the give box, optionally rotating the results to
# create a non-axis aligned distribution.
#===============================================================================
def randomPoints(numPoints,
                 xmin, xmax,
                 ymin, ymax,
                 theta = None):
    result = vector_of_Vector()

    # Determine the rotational transform.
    if theta is None:
        theta = rangen.uniform(0.0, 2.0*pi)
    R = rotationMatrix(Vector(cos(theta), sin(theta)))
    
    for i in xrange(numPoints):
        result.append(R*Vector(rangen.uniform(xmin, xmax),
                               rangen.uniform(ymin, ymax)))

    return R, result

#===============================================================================
# A local class for sorting points in a counter-clockwise fashion around a
# point.
#===============================================================================
class SortCounterClockwise:
    def __init__(self, p0):
        self.p0 = p0
        return
    def __call__(self, p1, p2):
        t = (p1 - self.p0).cross(p2 - self.p0).z
        if t > 0.0:
            return 1
        elif t < 0.0:
            return -1
        else:
            return 0

#===============================================================================
# Test class
#===============================================================================
class TestPolygon(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.ntests = 5000
        self.npoints = 1000
        self.xmin, self.xmax = -5.0, 2.0
        self.ymin, self.ymax = -1.0, 3.0
        self.R, self.points = randomPoints(self.npoints,
                                           self.xmin, self.xmax,
                                           self.ymin, self.ymax)
        self.polygon = Polygon(self.points)
        #self.plot = append(plotPolygon(self.polygon, False, True))
        return

    #---------------------------------------------------------------------------
    # Find the inner & outer radii of the polygon.
    #---------------------------------------------------------------------------
    def innerOuterRadii(self, polygon):
        rinner = 1e10
        router = 0.0
        centroid = polygon.centroid()
        for f in self.polygon.facets():
            r = abs((f.point1 - centroid).dot(f.normal))
            rinner = min(rinner, r)
        for v in polygon.vertices():
            r = (v - centroid).magnitude()
            router = max(router, r)
        router *= 1.0 + 1.0e-5
        return rinner, router

    #---------------------------------------------------------------------------
    # Check that all the seed points are contained in the polygon.
    #---------------------------------------------------------------------------
    def testContainSeeds(self):
        for p in self.points:
            self.failUnless(self.polygon.contains(p),
                            "Polygon does not contain seed point: %s" % str(p))
        return

    #---------------------------------------------------------------------------
    # Check that all the seed points are contained in the polygon using generic
    # contain method.
    #---------------------------------------------------------------------------
    def testGenericContainSeeds(self):
        for p in self.points:
            self.failUnless(pointInPolygon(p, self.polygon, True),
                            "Polygon does not contain seed point: %s" % str(p))
        return

    #---------------------------------------------------------------------------
    # Generate random points in the polygon and test they are contained.
    #---------------------------------------------------------------------------
    def testRandomInnerPoints(self):
        rinner, router = self.innerOuterRadii(self.polygon)
        centroid = self.polygon.centroid()
        for i in xrange(self.ntests):
            theta = rangen.uniform(0.0, 2.0*pi)
            p = centroid + rangen.uniform(0.0, rinner) * Vector(cos(theta), sin(theta))
            self.failUnless(self.polygon.contains(p),
                            "Polygon should contain %s but reports it does not." % str(p))
        return

    #---------------------------------------------------------------------------
    # Generate random points outside the polygon and test they are not contained.
    #---------------------------------------------------------------------------
    def testRandomOuterPoints(self):
        rinner, router = self.innerOuterRadii(self.polygon)
        centroid = self.polygon.centroid()
        for i in xrange(self.ntests):
            theta = rangen.uniform(0.0, 2.0*pi)
            p = centroid + rangen.uniform(router, 2.0*router) * Vector(cos(theta), sin(theta))
            self.failUnless(not self.polygon.contains(p),
                            "%s should be outside polygon but polygon reports it is contained." % str(p))
        return

    #---------------------------------------------------------------------------
    # Test vertex containment.
    #---------------------------------------------------------------------------
    def testVertexPointContainment(self):
        vertices = self.polygon.vertices()
        for v in vertices:
            self.failUnless(self.polygon.contains(v),
                            "%s vertex position should be contained." % str(v))
            self.failUnless(not self.polygon.contains(v, False),
                            "%s vertex position should not be contained." % str(v))
        return

    #---------------------------------------------------------------------------
    # Test that nested polygons intersect.
    #---------------------------------------------------------------------------
    def testIntersectInterior(self):
        centroid = self.polygon.centroid()
        vertices = vector_of_Vector(self.polygon.vertices())
        numVerts = len(vertices)
        for i in xrange(numVerts):
            vertices[i] = 0.5*(centroid + vertices[i])
            assert self.polygon.contains(vertices[i])
        polygon2 = Polygon(vertices)
        self.failUnless(self.polygon.intersect(polygon2),
                        "Failed to intersect with a contained polygon.")
        self.failUnless(polygon2.intersect(self.polygon),
                        "Failed to intersect when entirely contained within a polygon.")
        return
        
    #---------------------------------------------------------------------------
    # Test that two polygons just touching intersect.
    #---------------------------------------------------------------------------
    def testIntersectTouchingPolygons(self):
        xmin = self.polygon.xmin.x
        vertices = vector_of_Vector(self.polygon.vertices())
        numVerts = len(vertices)
        for i in xrange(numVerts):
            xv = vertices[i].x
            vertices[i].x -= 2.0*(xv - xmin)
        polygon2 = Polygon(vertices)
        self.failUnless(self.polygon.intersect(polygon2),
                        "Failed to intersect with polygon touching at a vertex")
        self.failUnless(polygon2.intersect(self.polygon),
                        "Failed to intersect with polygon touching at a vertex")
        return

    #---------------------------------------------------------------------------
    # Test that two separate polygons do not intersect.
    #---------------------------------------------------------------------------
    def testNotIntersectPolygons(self):
        xmin = self.polygon.xmin.x
        xlength = self.polygon.xmax.x - self.polygon.xmin.x
        vertices = vector_of_Vector(self.polygon.vertices())
        numVerts = len(vertices)
        for i in xrange(numVerts):
            xv = vertices[i].x
            vertices[i].x -= 2.0*(xv - xmin) + xlength
        polygon2 = Polygon(vertices)
        #plots.append(plotPolygon(self.polygon))
        #plots.append(plotPolygon(polygon2))
        self.failUnless(not self.polygon.intersect(polygon2),
                        "Erroneously claiming polygons intersect : [[%g,%g], [%g,%g], [[%g,%g], [%g,%g]]" % (self.polygon.xmin.x, 
                                                                                                             self.polygon.xmax.x,
                                                                                                             self.polygon.xmin.y,
                                                                                                             self.polygon.xmax.y,
                                                                                                             polygon2.xmin.x, 
                                                                                                             polygon2.xmax.x,
                                                                                                             polygon2.xmin.y,
                                                                                                             polygon2.xmax.y))
        self.failUnless(not polygon2.intersect(self.polygon),
                        "Erroneously claiming polygons intersect : [[%g,%g], [%g,%g], [[%g,%g], [%g,%g]]" % (self.polygon.xmin.x, 
                                                                                                             self.polygon.xmax.x,
                                                                                                             self.polygon.xmin.y,
                                                                                                             self.polygon.xmax.y,
                                                                                                             polygon2.xmin.x, 
                                                                                                             polygon2.xmax.x,
                                                                                                             polygon2.xmin.y,
                                                                                                             polygon2.xmax.y))
        return

    #---------------------------------------------------------------------------
    # Test that a nested box intersects.
    #---------------------------------------------------------------------------
    def testIntersectBoxInterior(self):
        rinner, router = self.innerOuterRadii(self.polygon)
        centroid = self.polygon.centroid()
        box = pair_Vector2d_Vector2d(centroid - 0.95*rinner*Vector.one,
                                     centroid + 0.95*rinner*Vector.one)
        self.failUnless(self.polygon.intersect(box),
                        "Failed to intersect with a contained box.")
        return
        
    #---------------------------------------------------------------------------
    # Test that a circumscribing box intersects.
    #---------------------------------------------------------------------------
    def testIntersectBoxExterior(self):
        rinner, router = self.innerOuterRadii(self.polygon)
        centroid = self.polygon.centroid()
        box = pair_Vector2d_Vector2d(centroid - 2.0*router*Vector.one,
                                     centroid + 2.0*router*Vector.one)
        self.failUnless(self.polygon.intersect(box),
                        "Failed to intersect with a box we are in.")
        return
        
    #---------------------------------------------------------------------------
    # Test that a box just touching a polygon intersects.
    #---------------------------------------------------------------------------
    def testIntersectTouchingBox(self):
        xmin = self.polygon.xmin.x
        vertex = None
        for v in self.polygon.vertices():
            if fuzzyEqual(v.x, xmin, 1.0e-10):
                vertex = v
        assert not vertex is None
        box = pair_Vector2d_Vector2d(Vector(v.x - 2.0, v.y - 2.0), v)
        self.failUnless(self.polygon.intersect(box),
                        "Failed to intersect with a box touching on one side.")

    #---------------------------------------------------------------------------
    # Test that that we don't intersect a box outside the polygon.
    #---------------------------------------------------------------------------
    def testNotIntersectBox(self):
        xmin, xmax = self.polygon.xmin, self.polygon.xmax
        delta = xmax - xmin
        box = pair_Vector2d_Vector2d(xmin - 2.0*delta, xmax - 2.0*delta)
        self.failUnless(not self.polygon.intersect(box),
                        "Erroneously intersect external box.")

    #---------------------------------------------------------------------------
    # Test reconstructing.
    #---------------------------------------------------------------------------
    def testReconstruct(self):
        polygon2 = Polygon()
        polygon2.reconstruct(self.polygon.vertices(),
                             self.polygon.facetVertices)
        self.failUnless(polygon2 == self.polygon,
                        "Failed to properly reconstruct polygon from vertices and facets.")
        return

    #---------------------------------------------------------------------------
    # Test volume
    #---------------------------------------------------------------------------
    def testVolume(self):
        verts0 = self.polygon.vertices()
        c = self.polygon.centroid()
        cmpmethod = SortCounterClockwise(c)
        verts = sorted(list(verts0), cmpmethod)
        p0 = verts[0]
        answer = 0.0
        for i in xrange(2, len(verts)):
            answer += (verts[i] - p0).cross(verts[i - 1] - p0).z
        answer *= 0.5
        vertstring = [str(x) for x in verts]
        self.failUnless(fuzzyEqual(self.polygon.volume, answer, 1.0e-10),
                        "Failed volume computation: %g != %g\n verts = %s" % (self.polygon.volume,
                                                                              answer,
                                                                              vertstring))

    #---------------------------------------------------------------------------
    # Closest point to vertices.
    #---------------------------------------------------------------------------
    def testClosestPointToVertices(self):
        verts = self.polygon.vertices()
        for p in verts:
            cp = self.polygon.closestPoint(p)
            self.failUnless(fuzzyEqual((cp - p).magnitude(), 0.0, 1.0e-10),
                            "Closest point to vertex %s : %s" % (p, cp))

    #---------------------------------------------------------------------------
    # Closest point on facets.
    #---------------------------------------------------------------------------
    def testClosestPointOnFacets(self):
        facets = self.polygon.facets()
        for f in facets:
            p = f.position
            cp = self.polygon.closestPoint(p)
            self.failUnless(fuzzyEqual((cp - p).magnitude(), 0.0, 1.0e-10),
                            "Closest point to facet position %s : %s" % (p, cp))

    #---------------------------------------------------------------------------
    # Closest point above facets.
    #---------------------------------------------------------------------------
    def testClosestPointAboveFacets(self):
        facets = self.polygon.facets()
        for f in facets:
            chi = rangen.uniform(0.1, 10.0)
            cp0 = f.position
            p = cp0 + chi*f.normal
            cp = self.polygon.closestPoint(p)
            self.failUnless(fuzzyEqual((cp0 - cp).magnitude(), 0.0, 1.0e-10),
                            "Closest point to position off of facet position %s : %s" % (cp0, cp))

    #---------------------------------------------------------------------------
    # Test ==
    #---------------------------------------------------------------------------
    def testEqual(self):
        self.failUnless(self.polygon == self.polygon,
                        "Failed self equivalence.")

    #---------------------------------------------------------------------------
    # Test !=
    #---------------------------------------------------------------------------
    def testNotEqual(self):
        polygon2 = Polygon()
        self.failUnless(polygon2 != self.polygon,
                        "Failed not equal.")

    #---------------------------------------------------------------------------
    # Test copy constructor
    #---------------------------------------------------------------------------
    def testCopy(self):
        polygon2 = Polygon(self.polygon)
        self.failUnless(polygon2 == self.polygon,
                        "Failed to copy construct.")

    #---------------------------------------------------------------------------
    # Test shift in-place
    #---------------------------------------------------------------------------
    def testShiftInPlace(self):
        shift = Vector(rangen.uniform(-10.0, -10.0),
                       rangen.uniform(-10.0, -10.0))
        polygon2 = Polygon(self.polygon)
        polygon2 += shift
        for p0, p1 in zip([self.polygon.xmin, self.polygon.xmax] + list(self.polygon.vertices()),
                          [polygon2.xmin, polygon2.xmax] + list(polygon2.vertices())):
            pshift = p0 + shift
            self.failUnless(pshift == p1, "In-place shift point comparison failed: %s != %s" % (pshift, p1))

if __name__ == "__main__":
    unittest.main()
