#ATS:test(SELF, label="Field unit tests")
from math import *
import unittest
from Spheral1d import *
from SpheralTestUtilities import uniqueSequence

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
        for i in range(n):
            result[i] = l[i]
    except:
        result = vector_of_int(l)    # pybind11
    return result

#-------------------------------------------------------------------------------
# Define our unit test class.
#-------------------------------------------------------------------------------
class testField(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def setUp(self):
        self.nodes = makeFluidNodeList("test bed",
                                       eos,
                                       n)
        self.field = IntField("test field", self.nodes)
        for i in range(n):
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
    # deleteElements (random)
    #---------------------------------------------------------------------------
    def testDeleteElementsRandom(self):
        elementsToKill = uniqueSequence([g.randint(0, n - 1) for i in range(n//10)])
        self.nodes.deleteNodes(vector_from_list(elementsToKill))
        answer = list(range(n))
        for i in elementsToKill:
            answer.remove(i)
        assert len(answer) == self.nodes.numInternalNodes
        for i in range(len(answer)):
            assert self.field[i] == answer[i]
        return

    #---------------------------------------------------------------------------
    # deleteElements (front)
    #---------------------------------------------------------------------------
    def testDeleteElementsFront(self):
        elementsToKill = list(range(n//10))
        self.nodes.deleteNodes(vector_from_list(elementsToKill))
        answer = list(range(n))
        for i in elementsToKill:
            answer.remove(i)
        assert len(answer) == self.nodes.numInternalNodes
        for i in range(len(answer)):
            assert self.field[i] == answer[i]
        return

    #---------------------------------------------------------------------------
    # deleteElements (back)
    #---------------------------------------------------------------------------
    def testDeleteElementsBack(self):
        elementsToKill = list(range(9*n//10, n))
        self.nodes.deleteNodes(vector_from_list(elementsToKill))
        answer = list(range(n))
        for i in elementsToKill:
            answer.remove(i)
        assert len(answer) == self.nodes.numInternalNodes
        for i in range(len(answer)):
            assert self.field[i] == answer[i]
        return

#-------------------------------------------------------------------------------
# Run the test
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
