#ATS:test(SELF, label="SiloFileIO unit tests")
from Spheral import *
from FileIOTestBase import *

import os
import unittest

#-------------------------------------------------------------------------------
# SiloFileIO tests.
#-------------------------------------------------------------------------------
class SiloFileIOTest(FileIOTestBase, unittest.TestCase):

    def setUp(self):
        self.n = 10 # 1000
        self.intmin = -2**24
        self.intmax = 2**24
        self.unsignedmin = 0
        self.unsignedmax = 2**32
        self.doublemin = -1e50
        self.doublemax = 1e50
        self.constructor = SiloFileIO

        # Size the NodeLists.
        nodes1d.numInternalNodes = self.n
        nodes2d.numInternalNodes = self.n
        nodes3d.numInternalNodes = self.n

        return

    def tearDown(self):
        return

    def removeFile(self, filename):
        os.remove(filename + ".silo")

    #---------------------------------------------------------------------------
    # Compoundarray
    #---------------------------------------------------------------------------
    def testCompoundarray(self):

        db = silo.DBCreate("TestCompoundarray.silo", 
                           silo.DB_CLOBBER, silo.DB_LOCAL, "some file", silo.DB_HDF5)
        thpt = vector_of_vector_of_int([vector_of_int(range(100)), vector_of_int(range(10))])
        elemNames = vector_of_string(["range(100)", "range(10)"])
        opts = silo.DBoptlist(1024)
        assert opts.addOption(silo.DBOPT_CYCLE, 10) == 0
        assert opts.addOption(silo.DBOPT_DTIME, 100.0) == 0
        assert silo.DBPutCompoundarray(db, "stuff", elemNames, thpt, opts) == 0
        assert silo.DBClose(db) == 0
        self.removeFile("TestCompoundarray")
        return

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    #print "waiting..."
    #raw_input()
    unittest.main()
