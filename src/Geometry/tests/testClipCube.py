#ATS:test(SELF, label="Polyhedron cube clipping tests")

import unittest
from math import *
import time
from SpheralTestUtilities import fuzzyEqual

from Spheral3d import *

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyhedronClipping(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.points = vector_of_Vector()
        for tup in [(0,0,0),
                    (1,0,0),
                    (0,1,0),
                    (1,1,0),
                    (0,0,1),
                    (1,0,1),
                    (0,1,1),
                    (1,1,1)]:
            self.points.append(Vector(*tup))
        self.cube = Polyhedron(self.points)
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing through the cube.
    #---------------------------------------------------------------------------
    def testClipInternal(self):
        planes1, planes2 = vector_of_Plane(), vector_of_Plane()
        planes1.append(Plane3d(Vector(0.5, 0.5, 0.5), Vector( 1,0,0)))
        planes2.append(Plane3d(Vector(0.5, 0.5, 0.5), Vector(-1,0,0)))
        chunk1 = Polyhedron(self.cube)
        chunk2 = Polyhedron(self.cube)
        from PolyhedronFileUtilities import writePolyhedronOBJ
        writePolyhedronOBJ(self.cube, "cube.obj")
        clipFacetedVolumeByPlanes(chunk1, planes1)
        clipFacetedVolumeByPlanes(chunk2, planes2)
        writePolyhedronOBJ(chunk1, "chunk_ONE.obj")
        writePolyhedronOBJ(chunk2, "chunk_TWO.obj")
        return

if __name__ == "__main__":
    unittest.main()
