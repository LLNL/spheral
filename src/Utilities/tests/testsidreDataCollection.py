from math import *
import unittest
from Spheral1d import *

import random
g = random.Random()

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
n = 1000
eos = GammaLawGasCGS(5.0/3.0, 2.0)
WT = TableKernel(BSplineKernel(), 1000)

#-------------------------------------------------------------------------------
#  Helper method to convert a list to vector_of_int
#-------------------------------------------------------------------------------
def vector_from_list(l):
    n = len(l)
    try:
        result = vector_of_int(n)    # pybindgen
        for i in xrange(n):
            result[i] = l[i]
    except:
        result = vector_of_int(l)    # pybind11
    return result

#-------------------------------------------------------------------------------
# Define our unit test class.
#-------------------------------------------------------------------------------
class testsidreDataCollection(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def setUp(self):
        self.nodes = makeFluidNodeList("test bed",
                                       eos,
                                       n)
        self.field = IntField("test field", self.nodes)
        for i in xrange(n):
            self.field[i] = i
        return

    #---------------------------------------------------------------------------
    # Destructor
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.field
        del self.nodes
        return

    #---------------------------------------------------------------------------
    # alloc_view (Field of int)
    #---------------------------------------------------------------------------
    def testAlloc_viewInt(self):
        alloc_view("IntSidreTest", self.nodes)
        return



#-------------------------------------------------------------------------------
# Run the test
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()