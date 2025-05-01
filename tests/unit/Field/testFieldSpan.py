#ATS:test(SELF, label="Field unit tests")
from math import *
import unittest
from Spheral2d import *
from SpheralTestUtilities import uniqueSequence

nInt = 1000
nGhost = 200
nTot = nInt + nGhost
nodes = makeVoidNodeList("test bed",
                         numInternal = nInt,
                         numGhost = nGhost)

#-------------------------------------------------------------------------------
# Set up some random generators
#-------------------------------------------------------------------------------
import random
g = random.Random()

def randInt():
    return g.randint(-10000, 100000)

def randScalar():
    return g.uniform(-1e10, 1e10)

def randTensor():
    return Tensor(randScalar(), randScalar(),
                  randScalar(), randScalar())

def randTRT():
    result = ThirdRankTensor()
    for j in range(len(result)):
        result[j] = randScalar()
    return result

#-------------------------------------------------------------------------------
# Test class
#-------------------------------------------------------------------------------
class testFieldSpan(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def setUp(self):
        self.ntests = 100

        self.intField = IntField("integer field", nodes)
        self.scalarField = ScalarField("scalar field", nodes)
        self.tensorField = TensorField("tensor field", nodes)
        self.trtField = ThirdRankTensorField("third rank tensor field", nodes)

        self.intSpan = IntFieldSpan(self.intField)
        self.scalarSpan = ScalarFieldSpan(self.scalarField)
        self.tensorSpan = TensorFieldSpan(self.tensorField)
        self.trtSpan = ThirdRankTensorFieldSpan(self.trtField)

        self.stuff = [(self.intField,    self.intSpan,      randInt),
                      (self.scalarField, self.scalarSpan,   randScalar),
                      (self.tensorField, self.tensorSpan,   randTensor),
                      (self.trtField,    self.trtSpan,      randTRT)]

        return

    #---------------------------------------------------------------------------
    # Destructor
    #---------------------------------------------------------------------------
    def tearDown(self):
        for field, span, generator in self.stuff:
            del field, span
        return

    #---------------------------------------------------------------------------
    # index (field -> span)
    #---------------------------------------------------------------------------
    def testIndexInternalField2Span(self):
        for k in range(self.ntests):
            i = g.randint(0, nTot - 1)
            for f, s, gen in self.stuff:
                f[i] = gen()
                assert f[i] == s[i]
        return

    #---------------------------------------------------------------------------
    # index (span -> field)
    #---------------------------------------------------------------------------
    def testIndexInternalSpan2Field(self):
        for k in range(self.ntests):
            i = g.randint(0, nTot - 1)
            for f, s, gen in self.stuff:
                s[i] = gen()
                assert f[i] == s[i]
        return

#-------------------------------------------------------------------------------
# Run the test
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
