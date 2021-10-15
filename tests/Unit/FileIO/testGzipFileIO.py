#ATS:test(SELF, label="GzipFileIO unit tests")
from Spheral import *
from FileIOTestBase import *

import os
import unittest

#-------------------------------------------------------------------------------
# GzipFileIO tests.
#-------------------------------------------------------------------------------
class GzipFileIOTest(FileIOTestBase, unittest.TestCase):

    def setUp(self):
        self.n = 10 # 1000
        self.intmin = -2**24
        self.intmax = 2**24
        self.unsignedmin = 0
        self.unsignedmax = 2**32
        self.doublemin = -1e50
        self.doublemax = 1e50
        self.constructor = GzipFileIO

        # Size the NodeLists.
        nodes1d.numInternalNodes = self.n
        nodes2d.numInternalNodes = self.n
        nodes3d.numInternalNodes = self.n

        return

    def tearDown(self):
        return

    def removeFile(self, filename):
        os.remove(filename + ".gz")

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
