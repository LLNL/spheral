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

# Make a cube
points = vector_of_Vector()
for coords in [(0,0,0), (1,0,0), (0,1,0), (1,1,0),
               (0,0,1), (1,0,1), (0,1,1), (1,1,1)]:
    points.append(Vector(*coords))
cube = Polyhedron(points)

# Make a non-convex notched thingy.
points = vector_of_Vector()
for coords in [(1,0,0), (1,4,0), (1,4,2), (1,3,2), (1,2,1), (1,1,2), (1,0,2),
               (0,0,0), (0,4,0), (0,4,2), (0,3,2), (0,2,1), (0,1,2), (0,0,2)]:
    points.append(Vector(*coords))
facets = vector_of_vector_of_unsigned()
for fac in [(0, 1, 2, 3, 4, 5, 6),
            (8, 7, 13, 12, 11, 10, 9),
            (0, 7, 8, 1),
            (1, 8, 9, 2),
            (0, 6, 13, 7),
            (6, 5, 12, 13),
            (5, 4, 11, 12),
            (4, 3, 10, 11),
            (3, 2, 9, 10)]:
    face = vector_of_unsigned()
    for i in fac:
        face.append(i)
    facets.append(face)
notchedthing = Polyhedron(points, facets)

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyhedronClipping(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.polyhedra = [cube, notchedthing]
        self.ntests = 1000
        return

    #---------------------------------------------------------------------------
    # Spheral::Polyhedron --> PolyClipper::Polyhedron
    #---------------------------------------------------------------------------
    def testConvertToPolyhedron(self):
        for poly in self.polyhedra:
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
        for poly0 in self.polyhedra:
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly0)
            assert poly0.vertices().size() == PCpoly.size()
            poly1 = Polyhedron()
            PolyClipper.convertFromPolyhedron(poly1, PCpoly)
            assert poly1.volume == poly1.volume
            assert poly1.centroid() == poly1.centroid()

    #---------------------------------------------------------------------------
    # PolyClipper::copyPolyhedron
    #---------------------------------------------------------------------------
    def testCopyPolyhedron(self):
        for poly0 in self.polyhedra:
            PCpoly0 = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly0, poly0)
            assert poly0.vertices().size() == PCpoly0.size()
            PCpoly1 = PolyClipper.Polyhedron()
            PolyClipper.copyPolyhedron(PCpoly1, PCpoly0)
            assert PCpoly0.size() == PCpoly1.size()
            vol0, centroid0 = PolyClipper.moments(PCpoly0)
            vol1, centroid1 = PolyClipper.moments(PCpoly1)
            assert vol0 == vol1
            assert centroid0 == centroid1

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polygon.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
        for poly in self.polyhedra:
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_Plane(), vector_of_Plane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane(p0,  phat))
                planes2.append(Plane(p0, -phat))
                PCchunk1 = PolyClipper.Polyhedron()
                PCchunk2 = PolyClipper.Polyhedron()
                PolyClipper.copyPolyhedron(PCchunk1, PCpoly)
                PolyClipper.copyPolyhedron(PCchunk2, PCpoly)
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
        for poly in self.polyhedra:
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                planes1, planes2 = vector_of_Plane(), vector_of_Plane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                phat = Vector(rangen.uniform(-1.0, 1.0), 
                              rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane(p0,  phat))
                planes2.append(Plane(p0,  phat))
                planes2.append(Plane(p0,  phat))
                PCchunk1 = PolyClipper.Polyhedron()
                PCchunk2 = PolyClipper.Polyhedron()
                PolyClipper.copyPolyhedron(PCchunk1, PCpoly)
                PolyClipper.copyPolyhedron(PCchunk2, PCpoly)
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
        for poly in self.polyhedra:
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes = vector_of_Plane()
                planes.append(Plane(p0, -phat))
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
        for poly in self.polyhedra:
            for i in xrange(self.ntests):
                planes = vector_of_Plane()
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector(cos(theta), sin(theta))
                p0 = poly.centroid() + r*phat
                planes.append(Plane(p0, phat))
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
        for poly in self.polyhedra:
            PCpoly = PolyClipper.Polyhedron()
            PolyClipper.convertToPolyhedron(PCpoly, poly)
            for i in xrange(self.ntests):
                planes1 = vector_of_Plane()
                p0 = Vector(rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0),
                            rangen.uniform(0.0, 1.0))
                for iplane in xrange(2):
                    planes1.append(Plane(point = p0,
                                         normal = Vector(rangen.uniform(-1.0, 1.0), 
                                                         rangen.uniform(-1.0, 1.0), 
                                                         rangen.uniform(-1.0, 1.0)).unitVector()))
                planes2 = vector_of_Plane(planes1)
                planes3 = vector_of_Plane(planes1)
                planes4 = vector_of_Plane(planes1)
                planes2[0].normal = -planes2[0].normal
                planes3[1].normal = -planes3[1].normal
                planes4[0].normal = -planes4[0].normal
                planes4[1].normal = -planes4[1].normal
                PCchunk1 = PolyClipper.Polyhedron()
                PCchunk2 = PolyClipper.Polyhedron()
                PCchunk3 = PolyClipper.Polyhedron()
                PCchunk4 = PolyClipper.Polyhedron()
                PolyClipper.copyPolyhedron(PCchunk1, PCpoly)
                PolyClipper.copyPolyhedron(PCchunk2, PCpoly)
                PolyClipper.copyPolyhedron(PCchunk3, PCpoly)
                PolyClipper.copyPolyhedron(PCchunk4, PCpoly)
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


if __name__ == "__main__":
    unittest.main()
