#ATS:test(SELF, label="PolyClipper 2D (polygon) tests")

import unittest
from math import *
import time
from PolyhedronFileUtilities import *

from Spheral2d import *
from SpheralTestUtilities import fuzzyEqual

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
# Make a square
#   3    2
#   |----|
#   |    |
#   |----|
#   0    1
square_points = vector_of_Vector()
for coords in [(0,0), (10,0), (10,10), (0,10)]:
    square_points.append(Vector(*coords))

#-------------------------------------------------------------------------------
# Make a non-convex notched thingy.
#    6           5      3          2
#    |------------\    /-----------|
#    |             \  /            |
#    |              \/             |
#    |               4             |
#    |                             |
#    |------------------------------
#    0                             1
notched_points = vector_of_Vector()
for coords in [(0,0), (4,0), (4,2), (3,2), (2,1), (1,2), (0,2)]:
    notched_points.append(Vector(*coords))

#-------------------------------------------------------------------------------
# A degenerate square
#   4    3
#  5|----|
#   |    |
#   |----|2
#   0    1
degenerate_square_points = vector_of_Vector()
for coords in [(0,0), (1,0), (1,0), (1,1), (0,1), (0,1)]:
    degenerate_square_points.append(Vector(*coords))

#-------------------------------------------------------------------------------
# Compute the vertex neighbors assuming an ordered ring of the given size.
#-------------------------------------------------------------------------------
def vertexNeighbors(points):
    n = len(points)
    neighbors = vector_of_vector_of_int(n, vector_of_int(2))
    for i in xrange(n):
        neighbors[i][0] = (i - 1) % n
        neighbors[i][1] = (i + 1) % n
    return neighbors

#-------------------------------------------------------------------------------
# Compute the facets assuming an ordered ring of the given size.
#-------------------------------------------------------------------------------
def facets(points):
    n = len(points)
    facets = vector_of_vector_of_unsigned(n, vector_of_unsigned(2))
    for i in xrange(n):
        facets[i][0] = i
        facets[i][1] = (i + 1) % n
    return facets

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyClipper2d(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.convexPointSets = [square_points, degenerate_square_points]
        self.nonconvexPointSets = [notched_points]
        self.pointSets = self.convexPointSets + self.nonconvexPointSets
        self.ntests = 10000

    #---------------------------------------------------------------------------
    # initializePolygon
    #---------------------------------------------------------------------------
    def test_initializePolygon(self):
        for points in self.pointSets:
            poly = Polygon(points, facets(points))
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            vol, centroid = PolyClipper.moments(PCpoly)
            self.failUnless(vol == poly.volume,
                            "Volume comparison failure: %g != %g" % (vol, poly.volume))
            self.failUnless(centroid == poly.centroid(),
                            "Centroid comparison failure: %s != %s" % (centroid, poly.centroid()))

    #---------------------------------------------------------------------------
    # collapseDegenerates
    #---------------------------------------------------------------------------
    def test_collapseDegenerates(self):
        PCpoly0 = PolyClipper.Polygon()
        PolyClipper.initializePolygon(PCpoly0, degenerate_square_points, vertexNeighbors(degenerate_square_points))
        assert PCpoly0.size() == len(degenerate_square_points)
        PCpoly1 = PolyClipper.Polygon(PCpoly0)
        PolyClipper.collapseDegenerates(PCpoly1, 1.0e-10)
        assert PCpoly1.size() == 4
        vol0, centroid0 = PolyClipper.moments(PCpoly0)
        vol1, centroid1 = PolyClipper.moments(PCpoly1)
        assert vol1 == vol0
        assert centroid1 == centroid0

    #---------------------------------------------------------------------------
    # Spheral::Polygon --> PolyClipper::Polygon
    #---------------------------------------------------------------------------
    def testConvertToPolygon(self):
        for points in self.pointSets:
            poly = Polygon(points, facets(points))
            PCpoly = PolyClipper.Polygon()
            PolyClipper.convertToPolygon(PCpoly, poly)
            assert poly.vertices().size() == PCpoly.size()
            vol, centroid = PolyClipper.moments(PCpoly)
            self.failUnless(vol == poly.volume,
                            "Volume comparison failure: %g != %g" % (vol, poly.volume))
            self.failUnless(centroid == poly.centroid(),
                            "Centroid comparison failure: %s != %s" % (centroid, poly.centroid()))


    #---------------------------------------------------------------------------
    # PolyClipper::Polygon --> Spheral::Polygon
    #---------------------------------------------------------------------------
    def testConvertFromPolygon(self):
        for points in self.pointSets:
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            poly0 = Polygon(points, facets(points))
            poly1 = Polygon()
            PolyClipper.convertFromPolygon(poly1, PCpoly)
            assert poly1 == poly0

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polygon.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
        for points in self.pointSets:
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            poly = Polygon(points, facets(points))
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_PolyClipperPlane(), vector_of_PolyClipperPlane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(PolyClipper.PolyClipperPlane2d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane2d(p0, -phat))
                PCchunk1 = PolyClipper.Polygon(PCpoly)
                PCchunk2 = PolyClipper.Polygon(PCpoly)
                PolyClipper.clipPolygon(PCchunk1, planes1)
                PolyClipper.clipPolygon(PCchunk2, planes2)
                chunk1 = Polygon()
                chunk2 = Polygon()
                PolyClipper.convertFromPolygon(chunk1, PCchunk1)
                PolyClipper.convertFromPolygon(chunk2, PCchunk2)
                success = fuzzyEqual(chunk1.volume + chunk2.volume, poly.volume)
                if not success:
                    print "Failed on pass ", i
                    print "Plane: ", p0, phat
                    print "Poly:\n", poly
                    print "Chunk 1:\n ", chunk1
                    print "Chunk 2:\n ", chunk2
                    vol1, cent1 = PolyClipper.moments(PCchunk1)
                    vol2, cent2 = PolyClipper.moments(PCchunk2)
                    print "Vol check: %g + %g = %g" % (vol1, vol2, vol1 + vol2)
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk1, "chunk_ONE.obj")
                    writePolyhedronOBJ(chunk2, "chunk_TWO.obj")
                self.failUnless(success,
                                "Plane clipping summing to wrong volumes: %s + %s != %s" % (chunk1.volume,
                                                                                            chunk2.volume,
                                                                                            poly.volume))

    #---------------------------------------------------------------------------
    # Clip with the same plane repeatedly.
    #---------------------------------------------------------------------------
    def testRedundantClip(self):
        for points in self.pointSets:
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            poly = Polygon(points, facets(points))
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_PolyClipperPlane(), vector_of_PolyClipperPlane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(PolyClipper.PolyClipperPlane2d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane2d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane2d(p0,  phat))
                PCchunk1 = PolyClipper.Polygon(PCpoly)
                PCchunk2 = PolyClipper.Polygon(PCpoly)
                PolyClipper.clipPolygon(PCchunk1, planes1)
                PolyClipper.clipPolygon(PCchunk2, planes2)
                chunk1 = Polygon()
                chunk2 = Polygon()
                PolyClipper.convertFromPolygon(chunk1, PCchunk1)
                PolyClipper.convertFromPolygon(chunk2, PCchunk2)
                success = fuzzyEqual(chunk1.volume, chunk2.volume)
                if not success:
                    print "Failed on pass ", i
                    print "Plane: ", p0, phat
                    print "Poly:\n", poly
                    print "Chunk 1:\n ", chunk1
                    print "Chunk 2:\n ", chunk2
                    vol1, cent1 = PolyClipper.moments(PCchunk1)
                    vol2, cent2 = PolyClipper.moments(PCchunk2)
                    print "Vol check: %g = %g" % (vol1, vol2)
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk1, "chunk_ONE.obj")
                    writePolyhedronOBJ(chunk2, "chunk_TWO.obj")
                self.failUnless(success,
                                "Redundant plane clipping wrong volumes: %s != %s" % (chunk1.volume,
                                                                                      chunk2.volume))

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polygon -- null test.
    #---------------------------------------------------------------------------
    def testNullClipOnePlane(self):
        for points in self.pointSets:
            poly = Polygon(points, facets(points))
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes = vector_of_PolyClipperPlane()
                planes.append(PolyClipper.PolyClipperPlane2d(p0, -phat))
                PCchunk = PolyClipper.Polygon()
                PolyClipper.convertToPolygon(PCchunk, poly)
                PolyClipper.clipPolygon(PCchunk, planes)
                chunk = Polygon()
                PolyClipper.convertFromPolygon(chunk, PCchunk)
                success = (chunk.volume == poly.volume)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk, "chunk.obj")
                self.failUnless(success,
                                "Null plane clipping failure: %s != %s" % (chunk.volume, poly.volume))

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polygon and rejecting the whole thing.
    #---------------------------------------------------------------------------
    def testFullClipOnePlane(self):
        for points in self.pointSets:
            poly = Polygon(points, facets(points))
            for i in xrange(self.ntests):
                planes = vector_of_PolyClipperPlane()
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes.append(PolyClipper.PolyClipperPlane2d(p0, phat))
                PCchunk = PolyClipper.Polygon()
                PolyClipper.convertToPolygon(PCchunk, poly)
                PolyClipper.clipPolygon(PCchunk, planes)
                chunk = Polygon()
                PolyClipper.convertFromPolygon(chunk, PCchunk)
                success = (chunk.volume == 0.0)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk, "chunk.obj")
                self.failUnless(success,
                                "Full plane clipping failure: %s != %s" % (chunk.volume, poly.volume))

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polygon.
    #---------------------------------------------------------------------------
    def testClipInternalTwoPlanes(self):
        for points in self.pointSets:
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            poly = Polygon(points, facets(points))
            for i in xrange(self.ntests):
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                norm1 = Vector(rangen.uniform(-1.0, 1.0), 
                               rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector(rangen.uniform(-1.0, 1.0), 
                               rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = vector_of_PolyClipperPlane()
                planes1.append(PolyClipper.PolyClipperPlane2d(p0,  norm1))
                planes1.append(PolyClipper.PolyClipperPlane2d(p0,  norm2))
                planes2 = vector_of_PolyClipperPlane()
                planes2.append(PolyClipper.PolyClipperPlane2d(p0,  norm1))
                planes2.append(PolyClipper.PolyClipperPlane2d(p0, -norm2))
                planes3 = vector_of_PolyClipperPlane()
                planes3.append(PolyClipper.PolyClipperPlane2d(p0, -norm1))
                planes3.append(PolyClipper.PolyClipperPlane2d(p0,  norm2))
                planes4 = vector_of_PolyClipperPlane()
                planes4.append(PolyClipper.PolyClipperPlane2d(p0, -norm1))
                planes4.append(PolyClipper.PolyClipperPlane2d(p0, -norm2))
                PCchunk1 = PolyClipper.Polygon(PCpoly)
                PCchunk2 = PolyClipper.Polygon(PCpoly)
                PCchunk3 = PolyClipper.Polygon(PCpoly)
                PCchunk4 = PolyClipper.Polygon(PCpoly)
                PolyClipper.clipPolygon(PCchunk1, planes1)
                PolyClipper.clipPolygon(PCchunk2, planes2)
                PolyClipper.clipPolygon(PCchunk3, planes3)
                PolyClipper.clipPolygon(PCchunk4, planes4)
                chunk1 = Polygon(poly)
                chunk2 = Polygon(poly)
                chunk3 = Polygon(poly)
                chunk4 = Polygon(poly)
                PolyClipper.convertFromPolygon(chunk1, PCchunk1)
                PolyClipper.convertFromPolygon(chunk2, PCchunk2)
                PolyClipper.convertFromPolygon(chunk3, PCchunk3)
                PolyClipper.convertFromPolygon(chunk4, PCchunk4)
                success = fuzzyEqual(chunk1.volume + chunk2.volume + chunk3.volume + chunk4.volume, poly.volume)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk1, "chunk_1ONE_TWOPLANES.obj")
                    writePolyhedronOBJ(chunk2, "chunk_2TWO_TWOPLANES.obj")
                    writePolyhedronOBJ(chunk3, "chunk_3THREE_TWOPLANES.obj")
                    writePolyhedronOBJ(chunk4, "chunk_4FOUR_TWOPLANES.obj")
                self.failUnless(success,
                                "Two plane clipping summing to wrong volumes: %s + %s + %s + %s = %s != %s" % (chunk1.volume,
                                                                                                               chunk2.volume,
                                                                                                               chunk3.volume,
                                                                                                               chunk4.volume,
                                                                                                               chunk1.volume + chunk2.volume + chunk3.volume + chunk4.volume,
                                                                                                               poly.volume))

    #---------------------------------------------------------------------------
    # Split a (convex) polygon into triangles.
    #---------------------------------------------------------------------------
    def testSplitIntoTriangles(self):
        for points in self.convexPointSets:
            PCpoly = PolyClipper.Polygon()
            PolyClipper.initializePolygon(PCpoly, points, vertexNeighbors(points))
            tris = PolyClipper.splitIntoTriangles(PCpoly)
            vol0, centroid0 = PolyClipper.moments(PCpoly)
            volTris = 0.0
            centroidTris = Vector()
            for inds in tris:
                assert len(inds) == 3
                a = ((PCpoly[inds[1]].position - PCpoly[inds[0]].position).cross(PCpoly[inds[2]].position - PCpoly[inds[0]].position).z)
                assert a >= 0.0
                volTris += a
                centroidTris += a*(PCpoly[inds[0]].position + PCpoly[inds[1]].position + PCpoly[inds[2]].position)
            volTris *= 0.5
            centroidTris /= 6.0*volTris
            assert abs(volTris - vol0) < 1.0e-20
            assert (centroidTris - centroid0).magnitude() < 1.0e-20

if __name__ == "__main__":
    unittest.main()
