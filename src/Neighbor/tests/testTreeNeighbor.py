#ATS:test(SELF, np=1, label="TreeNeighbor unit tests")
from math import *
import unittest
import random

from Spheral import *
from NeighborTestBase import *

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
class TestTreeNeighborRandom1d(NeighborRandom1d):
    NeighborRandom1d._NeighborType = TreeNeighbor1d

# #===============================================================================
# # Radom node distribution -- 2-D.
# #===============================================================================
# class TestTreeNeighborRandom2d(NeighborRandom2d):
#     NeighborRandom2d._NeighborType = TreeNeighbor2d

# #===============================================================================
# # Radom node distribution -- 3-D.
# #===============================================================================
# class TestTreeNeighborRandom3d(NeighborRandom3d):
#     NeighborRandom3d._NeighborType = TreeNeighbor3d

# #===============================================================================
# # Cylindrical node distribution -- 2-D.
# #===============================================================================
# class TestTreeNeighborCylindrical2d(NeighborCylindrical2d):
#     NeighborCylindrical2d._NeighborType = TreeNeighbor2d

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
