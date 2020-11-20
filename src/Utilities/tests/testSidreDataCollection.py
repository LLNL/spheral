from math import *
import unittest
from Spheral1d import *
from SpheralTestUtilities import *

# import sys

import string

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
n = 1000
eos = GammaLawGasCGS(5.0/3.0, 2.0)

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
class testSidreDataCollection(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def setUp(self):
        self.nodes = makeFluidNodeList("test bed",
                                       eos,
                                       n)
        # self.field = IntField("test field", self.nodes)
        self.SidreDataCollection = SidreDataCollection()
        # for i in xrange(n):
        #     self.field[i] = i
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
        self.field = IntField("test field", self.nodes)
        for i in xrange(n):
            self.field[i] = i
        
        answer = self.SidreDataCollection.alloc_view("IntSidreTest", self.field).getDataA(n)
        # for i in self.field:
        #     print(self.field[i]),
        # print(id(self.field))
        # print '[%s]' % ', '.join(map(str, answer))
        # print(answer[0])
        # for i in xrange(n):
        #     print(answer[i]),
        # self.SidreDataCollection.printDataStore()
        # assert sys.getsizeof(self.field[0]) == sys.getsizeof(answer[0])
        assert type(self.field[0]) == type(answer[0])
        assert len(self.field) == len(answer)
        for i in xrange(n):
            assert self.field[i] == answer[i]
        return
    
    # #---------------------------------------------------------------------------
    # # alloc_view (Field of char)
    # #---------------------------------------------------------------------------
    # def testAlloc_viewChar(self):
    #     self.field = ScalarField("test field", self.nodes)
    #     # for i in xrange(n):
    #     #     self.field[i] = string.ascii_letters[i%52]
    #     # for i in xrange(n):
    #     #     print(string.ascii_letters[i%52]),
    #     # print(string.ascii_letters[3])
    #     # print(chr(98))

        
    #     # answer = self.SidreDataCollection.alloc_view("CharSidreTest", self.field).getDataA(n)
    #     # assert sys.getsizeof(self.field[0]) == sys.getsizeof(answer[0])
    #     # assert len(self.field) == len(answer)
    #     # for i in xrange(n):
    #     #     assert self.field[i] == answer[i]
    #     return

    #---------------------------------------------------------------------------
    # alloc_view (Field of double)
    #---------------------------------------------------------------------------
    def testAlloc_viewDouble(self):
        self.field = ScalarField("double field", self.nodes)
        for i in xrange(n):
            self.field[i] = i

        print(type(self.field[0]))
        answer = self.SidreDataCollection.alloc_view("DoubleSidreTest", self.field).getDataB(n)
        # self.SidreDataCollection.printDataStore()
        # assert sys.getsizeof(self.field[0]) == sys.getsizeof(answer[0])
        assert type(self.field[0]) == type(answer[0])
        assert len(self.field) == len(answer)
        for i in xrange(n):
            assert self.field[i] == answer[i]
        return

    #---------------------------------------------------------------------------
    # alloc_view (Field of uint64_t)
    #---------------------------------------------------------------------------
    def testAlloc_viewUint64(self):
        self.field = ULLField("unsigned field", self.nodes)
        for i in xrange(n):
            self.field[i] = i
        print(type(self.field[0]))
        answer = self.SidreDataCollection.alloc_view("UnsignedSidreTest", self.field).getDataC(n)
        # self.SidreDataCollection.printDataStore()
        # assert sys.getsizeof(self.field[0]) == sys.getsizeof(answer[0])
        assert type(self.field[0]) == type(answer[0])
        assert len(self.field) == len(answer)
        for i in xrange(n):
            assert self.field[i] == answer[i]
        return

    #---------------------------------------------------------------------------
    # alloc_view (Field of vector of int)
    #---------------------------------------------------------------------------
    # def testAlloc_viewVectorInt(self):
    #     self.field = VectorIntField("test Vector Int field", self.nodes)
    #     for i in xrange(n):
    #         self.field[i] = vector_from_list(xrange(n))
        
    #     # answer = self.SidreDataCollection.alloc_view("IntSidreTest", self.field).getDataA(n)
    #     # for i in self.field:
    #     #     print(self.field[i]),
    #     # print(id(self.field))
    #     # print '[%s]' % ', '.join(map(str, answer))
    #     # print(answer[0])
    #     # for i in xrange(n):
    #     #     print(answer[i]),
    #     # self.SidreDataCollection.printDataStore()
    #     # assert sys.getsizeof(self.field[0]) == sys.getsizeof(answer[0])
    #     # assert len(self.field) == len(answer)
    #     # for i in xrange(n):
    #     #     assert self.field[i] == answer[i]
    #     return


#-------------------------------------------------------------------------------
# Run the test
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()