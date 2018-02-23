#ATS:test(SELF, label="Polyhedron clipping tests")

import unittest
from math import *
import time

from Spheral3d import *
from SpheralTestUtilities import fuzzyEqual
from PolyhedronFileUtilities import *

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
# Make a cube                   |y     
#                               |      
#                               |____x 
#   3/-----/2                  /       
#   /     /|                  /z       
# 7|-----|6|
#  |     | |
#  |  0  | /1
#  |_____|/
#  4     5
#
cube_points = vector_of_Vector()
for coords in [(0,0,0),  (10,0,0),  (10,10,0),  (0,10,0),
               (0,0,10), (10,0,10), (10,10,10), (0,10,10)]:
    cube_points.append(Vector(*coords))
cube_neighbors = vector_of_vector_of_int(8, vector_of_int(3))
for i, nghbs in enumerate([(1, 4, 3),
                           (5, 0, 2),
                           (3, 6, 1),
                           (7, 2, 0),
                           (5, 7, 0),
                           (1, 6, 4),
                           (5, 2, 7),
                           (4, 6, 3)]):
    for j, val in enumerate(nghbs):
        cube_neighbors[i][j] = val

cube_facets = vector_of_vector_of_unsigned(6, vector_of_unsigned(4))
for i, vinds in enumerate([(4, 5, 6, 7),
                           (1, 2, 6, 5),
                           (0, 3, 2, 1),
                           (4, 7, 3, 0),
                           (6, 2, 3, 7),
                           (1, 5, 4, 0)]):
    for j, val in enumerate(vinds):
        cube_facets[i][j] = val

#-------------------------------------------------------------------------------
# Make a non-convex notched thingy.                            |y     
#                                                              |      
#                                                              |____x 
#       6            5       3          2                     /       
#       /------------/      /----------/                     /z       
#      /            / \    /          /|
#     /            /   \ 4/          / |
#    /            /     \/          /  |
#    |------------\      /---------|9  |
#    |13        12 \    / 10       |   |
#    |              \  /           |   |
#    |               \/            |   |
#    |               11            |   |
#    |   0                         |  / 1
#    |                             | /
#    |------------------------------/
#    7                             8

notched_points = vector_of_Vector()
for coords in [(0,0,0), (4,0,0), (4,2,0), (3,2,0), (2,1,0), (1,2,0), (0,2,0),
               (0,0,1), (4,0,1), (4,2,1), (3,2,1), (2,1,1), (1,2,1), (0,2,1)]:
    notched_points.append(Vector(*coords))
notched_neighbors = vector_of_vector_of_int(14, vector_of_int(3))
for i, nghbs in enumerate([(7, 6, 1),   # 0
                           (0, 2, 8),   # 1
                           (1, 3, 9),   # 2
                           (4, 10, 2),  # 3
                           (5, 11, 3),  # 4
                           (6, 12, 4),  # 5
                           (13, 5, 0),  # 6
                           (8, 13, 0),  # 7
                           (1, 9, 7),   # 8
                           (2, 10, 8),  # 9
                           (9, 3, 11),  # 10
                           (10, 4, 12), # 11
                           (11, 5, 13), # 12
                           (7, 12, 6)]):# 13
    for j, val in enumerate(nghbs):
        notched_neighbors[i][j] = val

notched_facets = vector_of_vector_of_unsigned()
for fac in [(6, 5, 4, 3, 2, 1, 0),
            (7, 8, 9, 10, 11, 12, 13),
            (1, 2, 9, 8),
            (2, 3, 10, 9),
            (3, 4, 11, 10),
            (4, 5, 12, 11),
            (5, 6, 13, 12),
            (7, 13, 6, 0),
            (0, 1, 8, 7)]:
    face = vector_of_unsigned()
    for i in fac:
        face.append(i)
    notched_facets.append(face)

#-------------------------------------------------------------------------------
# A degenerate pyramid.  Just reuse the cube, but collapse one face.
degenerate_cube_points1 = vector_of_Vector()
for coords in [(0,0,0), (1,0,0), (1,1,0), (0,1,0),
               (0,0,1), (0,0,1), (0,0,1), (0,0,1)]:
    degenerate_cube_points1.append(Vector(*coords))

# Another one collapsing a different vertex.
degenerate_cube_points2 = vector_of_Vector()
for coords in [(0,0,0),  (10,0,0),  (10,10,0),  (10,10,0),
               (0,0,10), (10,0,10), (10,10,0), (10,10,0)]:
    degenerate_cube_points2.append(Vector(*coords))

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyhedronClipping(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.convexPolyData = [(cube_points, cube_neighbors, cube_facets),
                               (degenerate_cube_points1, cube_neighbors, cube_facets),
                               (degenerate_cube_points2, cube_neighbors, cube_facets)]
        self.nonconvexPolyData = [(notched_points, notched_neighbors, notched_facets)]
        self.degeneratePolyData = [(degenerate_cube_points1, cube_neighbors, cube_facets),
                                   (degenerate_cube_points2, cube_neighbors, cube_facets)]
        self.polyData = self.convexPolyData + self.nonconvexPolyData
        self.ntests = 1000
        return

    #---------------------------------------------------------------------------
    # initializePolyhedron
    #---------------------------------------------------------------------------
    def test_initializePolyhedron(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.initializePolyhedron(PCpoly, points, neighbors)
            vol, centroid = PolyClipper.moments(PCpoly)
            self.failUnless(vol == poly.volume,
                            "Volume comparison failure: %g != %g" % (vol, poly.volume))
            self.failUnless(centroid == poly.centroid(),
                            "Centroid comparison failure: %s != %s" % (centroid, poly.centroid()))

    #---------------------------------------------------------------------------
    # collapseDegenerates
    #---------------------------------------------------------------------------
    def test_collapseDegenerates(self):
        for points, neighbors, facets in self.degeneratePolyData:
            PCpoly0 = PolyClipper.Polyhedron()
            PolyClipper.initializePolyhedron(PCpoly0, points, neighbors)
            assert PCpoly0.size() == len(points)
            PCpoly1 = PolyClipper.Polyhedron(PCpoly0)
            PolyClipper.collapseDegenerates(PCpoly1, 1.0e-10)
            assert PCpoly1.size() == 5
            vol0, centroid0 = PolyClipper.moments(PCpoly0)
            vol1, centroid1 = PolyClipper.moments(PCpoly1)
            assert vol1 == vol0
            assert centroid1 == centroid0

    #---------------------------------------------------------------------------
    # Spheral::Polyhedron --> PolyClipper::Polyhedron
    #---------------------------------------------------------------------------
    def testConvertToPolyhedron(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            assert poly.vertices().size() == PCpoly.size()
            vol, centroid = PolyClipper.moments(PCpoly)
            self.failUnless(vol == poly.volume,
                            "Volume comparison failure: %g != %g" % (vol, poly.volume))
            self.failUnless(centroid == poly.centroid(),
                            "Centroid comparison failure: %s != %s" % (centroid, poly.centroid()))


    #---------------------------------------------------------------------------
    # PolyClipper::Polyhedron --> Spheral::Polyhedron
    #---------------------------------------------------------------------------
    def testConvertFromPolyhedron(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            assert poly.vertices().size() == PCpoly.size()
            poly1 = Polyhedron()
            PolyClipper.convertFromPolyhedron(poly1, PCpoly)
            assert poly1.volume == poly.volume
            assert poly1.centroid() == poly.centroid()

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polyhedron.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_PolyClipperPlane(), vector_of_PolyClipperPlane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(PolyClipper.PolyClipperPlane3d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane3d(p0, -phat))
                PCchunk1 = PolyClipper.Polyhedron(PCpoly)
                PCchunk2 = PolyClipper.Polyhedron(PCpoly)
                PolyClipper.clipPolyhedron(PCchunk1, planes1)
                PolyClipper.clipPolyhedron(PCchunk2, planes2)
                chunk1 = Polyhedron()
                chunk2 = Polyhedron()
                PolyClipper.convertFromPolyhedron(chunk1, PCchunk1)
                PolyClipper.convertFromPolyhedron(chunk2, PCchunk2)
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
        return

    #---------------------------------------------------------------------------
    # Clip with the same plane repeatedly.
    #---------------------------------------------------------------------------
    def testRedundantClip(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_PolyClipperPlane(), vector_of_PolyClipperPlane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(PolyClipper.PolyClipperPlane3d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane3d(p0,  phat))
                planes2.append(PolyClipper.PolyClipperPlane3d(p0,  phat))
                PCchunk1 = PolyClipper.Polyhedron(PCpoly)
                PCchunk2 = PolyClipper.Polyhedron(PCpoly)
                PolyClipper.clipPolyhedron(PCchunk1, planes1)
                PolyClipper.clipPolyhedron(PCchunk2, planes2)
                chunk1 = Polyhedron()
                chunk2 = Polyhedron()
                PolyClipper.convertFromPolyhedron(chunk1, PCchunk1)
                PolyClipper.convertFromPolyhedron(chunk2, PCchunk2)
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
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polyhedron -- null test.
    #---------------------------------------------------------------------------
    def testNullClipOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes = vector_of_PolyClipperPlane()
                planes.append(PolyClipper.PolyClipperPlane3d(p0, -phat))
                PCchunk = PolyClipper.Polyhedron()
                PolyClipper.convertToPolyhedron(PCchunk, poly)
                PolyClipper.clipPolyhedron(PCchunk, planes)
                chunk = Polyhedron()
                PolyClipper.convertFromPolyhedron(chunk, PCchunk)
                success = (chunk.volume == poly.volume)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk, "chunk.obj")
                self.failUnless(success,
                                "Null plane clipping failure: %s != %s" % (chunk.volume, poly.volume))
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polyhedron and rejecting the whole thing.
    #---------------------------------------------------------------------------
    def testFullClipOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            for i in xrange(self.ntests):
                planes = vector_of_PolyClipperPlane()
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes.append(PolyClipper.PolyClipperPlane3d(p0, phat))
                PCchunk = PolyClipper.Polyhedron()
                PolyClipper.convertToPolyhedron(PCchunk, poly)
                PolyClipper.clipPolyhedron(PCchunk, planes)
                chunk = Polyhedron()
                PolyClipper.convertFromPolyhedron(chunk, PCchunk)
                success = (chunk.volume == 0.0)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk, "chunk.obj")
                self.failUnless(success,
                                "Full plane clipping failure: %s != %s" % (chunk.volume, poly.volume))
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polyhedron.
    #---------------------------------------------------------------------------
    def testClipInternalTwoPlanes(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron(points, facets)
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                norm1 = Vector(rangen.uniform(-1.0, 1.0), 
                               rangen.uniform(-1.0, 1.0),
                               rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector(rangen.uniform(-1.0, 1.0), 
                               rangen.uniform(-1.0, 1.0),
                               rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = vector_of_PolyClipperPlane()
                planes1.append(PolyClipper.PolyClipperPlane3d(p0,  norm1))
                planes1.append(PolyClipper.PolyClipperPlane3d(p0,  norm2))
                planes2 = vector_of_PolyClipperPlane()
                planes2.append(PolyClipper.PolyClipperPlane3d(p0,  norm1))
                planes2.append(PolyClipper.PolyClipperPlane3d(p0, -norm2))
                planes3 = vector_of_PolyClipperPlane()
                planes3.append(PolyClipper.PolyClipperPlane3d(p0, -norm1))
                planes3.append(PolyClipper.PolyClipperPlane3d(p0,  norm2))
                planes4 = vector_of_PolyClipperPlane()
                planes4.append(PolyClipper.PolyClipperPlane3d(p0, -norm1))
                planes4.append(PolyClipper.PolyClipperPlane3d(p0, -norm2))
                PCchunk1 = PolyClipper.Polyhedron(PCpoly)
                PCchunk2 = PolyClipper.Polyhedron(PCpoly)
                PCchunk3 = PolyClipper.Polyhedron(PCpoly)
                PCchunk4 = PolyClipper.Polyhedron(PCpoly)
                PolyClipper.clipPolyhedron(PCchunk1, planes1)
                PolyClipper.clipPolyhedron(PCchunk2, planes2)
                PolyClipper.clipPolyhedron(PCchunk3, planes3)
                PolyClipper.clipPolyhedron(PCchunk4, planes4)
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                chunk3 = Polyhedron(poly)
                chunk4 = Polyhedron(poly)
                PolyClipper.convertFromPolyhedron(chunk1, PCchunk1)
                PolyClipper.convertFromPolyhedron(chunk2, PCchunk2)
                PolyClipper.convertFromPolyhedron(chunk3, PCchunk3)
                PolyClipper.convertFromPolyhedron(chunk4, PCchunk4)
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
        return

    #---------------------------------------------------------------------------
    # Split a (convex) polyhedron into tetrahedra.
    #---------------------------------------------------------------------------
    def testSplitIntoTetrahedra(self):
        for points, neighbors, facets in self.convexPolyData:
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.initializePolyhedron(PCpoly, points, neighbors)
            tets = PolyClipper.splitIntoTetrahedra(PCpoly)
            vol0, centroid0 = PolyClipper.moments(PCpoly)
            volTets = 0.0
            centroidTets = Vector()
            for inds in tets:
                assert len(inds) == 4
                v0 = PCpoly[inds[0]].position
                v1 = PCpoly[inds[1]].position
                v2 = PCpoly[inds[2]].position
                v3 = PCpoly[inds[3]].position
                V = (v1 - v0).dot((v2 - v0).cross(v3 - v0))
                assert V >= 0.0
                volTets += V
                centroidTets += V*(v0 + v1 + v2 + v3)
            volTets /= 6.0
            centroidTets /= 24.0*volTets
            assert abs(volTets - vol0) < 1.0e-20
            assert (centroidTets - centroid0).magnitude() < 1.0e-20

if __name__ == "__main__":
    unittest.main()
