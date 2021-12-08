#ATS:test(SELF, np=1, level=100, label="NestedGridNeighbor unit tests")
from math import *
import unittest
import random

from Spheral import *
from NeighborTestBase import *

# NestedGridNeighbor doesn't do ghost->ghost connectivity, so the overlap
# neighbor tests will choke.
del NeighborTestBase.testConnectivityMapOverlapNeighbors

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
NeighborRandom1d._NeighborType = NestedGridNeighbor1d

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
NeighborRandom2d._NeighborType = NestedGridNeighbor2d

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
NeighborRandom3d._NeighborType = NestedGridNeighbor3d

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
NeighborCylindrical2d._NeighborType = NestedGridNeighbor2d

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
