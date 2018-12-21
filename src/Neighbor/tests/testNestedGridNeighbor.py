#ATS:test(SELF, np=1, label="NestedGridNeighbor unit tests")
from math import *
import unittest
import random

from Spheral import *
from NeighborTestBase import *

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
class TestNestedGridNeighborRandom1d(NeighborRandom1d):
    NeighborRandom1d._NeighborType = NestedGridNeighbor1d

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
class TestNestedGridNeighborRandom2d(NeighborRandom2d):
    NeighborRandom2d._NeighborType = NestedGridNeighbor2d

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
class TestNestedGridNeighborRandom3d(NeighborRandom3d):
    NeighborRandom3d._NeighborType = NestedGridNeighbor3d

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
class TestNestedGridNeighborCylindrical2d(NeighborCylindrical2d):
    NeighborCylindrical2d._NeighborType = NestedGridNeighbor2d

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
