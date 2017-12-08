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
    # Clip with planes passing through the polyhedron.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
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
            for poly in self.polyhedra:
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                clipFacetedVolumeByPlanes(chunk1, planes1)
                clipFacetedVolumeByPlanes(chunk2, planes2)
                success = fuzzyEqual(chunk1.volume + chunk2.volume, poly.volume)
                if not success:
                    writePolyhedronOBJ(poly, "poly.obj")
                    writePolyhedronOBJ(chunk1, "chunk_ONE.obj")
                    writePolyhedronOBJ(chunk2, "chunk_TWO.obj")
                self.failUnless(success,
                                "Plane clipping summing to wrong volumes: %s + %s != %s" % (chunk1.volume,
                                                                                            chunk2.volume,
                                                                                            poly.volume))
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polyhedron -- null test.
    #---------------------------------------------------------------------------
    def testNullClipOnePlane(self):
        for poly in self.polyhedra:
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * (poly.xmax - poly.xmin).magnitude()
                theta = rangen.uniform(0.0, 2.0*pi)
                phi = rangen.uniform(0.0, pi)
                phat = Vector(cos(theta)*sin(phi),
                              sin(theta)*sin(phi),
                              sin(phi)).unitVector()
                p0 = poly.centroid() + r*phat
                planes = vector_of_Plane()
                planes.append(Plane(p0, -phat))
                chunk = Polyhedron(poly)
                clipFacetedVolumeByPlanes(chunk, planes)
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
                phi = rangen.uniform(0.0, pi)
                phat = Vector(cos(theta)*sin(phi),
                              sin(theta)*sin(phi),
                              sin(phi)).unitVector()
                p0 = poly.centroid() + r*phat
                planes.append(Plane(p0, phat))
                chunk = Polyhedron(poly)
                clipFacetedVolumeByPlanes(chunk, planes)
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
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                chunk3 = Polyhedron(poly)
                chunk4 = Polyhedron(poly)
                clipFacetedVolumeByPlanes(chunk1, planes1)
                clipFacetedVolumeByPlanes(chunk2, planes2)
                clipFacetedVolumeByPlanes(chunk3, planes3)
                clipFacetedVolumeByPlanes(chunk4, planes4)
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
