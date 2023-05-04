#ATS:test(SELF, label="SidreFileIO unit tests")
from Spheral import *
from FileIOTestBase import *

import os, shutil
import unittest
import mpi
import shutil

#-------------------------------------------------------------------------------
# SidreFileIO tests.
#-------------------------------------------------------------------------------
class SidreFileIOTest(FileIOTestBase, unittest.TestCase):

    def setUp(self):
        self.n = 10 # 1000
        self.intmin = -2**24
        self.intmax = 2**24
        self.unsignedmin = 0
        self.unsignedmax = 2**32
        self.doublemin = -1e10
        self.doublemax = 1e10
        self.constructor = SidreFileIO

        # Size the NodeLists.
        nodes1d.numInternalNodes = self.n
        nodes2d.numInternalNodes = self.n
        nodes3d.numInternalNodes = self.n

        return

    def tearDown(self):
        return

    # If we are using MPI then we need to remove a directory because we are using Spio,
    # otherwise we remove a file as is the case with the other FileIO types.
    def removeFile(self, filename):
        os.remove(filename)
        if not mpi.is_fake_mpi():
            shutil.rmtree(filename)

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    #print "waiting..."
    #raw_input()
    unittest.main()
