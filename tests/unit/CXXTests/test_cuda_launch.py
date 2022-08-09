#ATS:test(SELF, ngpu=1, label="Spheral CUDA Launch test.") 

from SpheralCXXTests import *
import unittest

class TestCUDALaunch(unittest.TestCase):

  def testCUDALaunch(self):
    self.assertEqual(launchCaller(3,4), 7, "Should be 7")
    return

if __name__ == "__main__":
  unittest.main()
