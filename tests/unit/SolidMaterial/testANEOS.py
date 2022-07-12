#ATS:test(SELF, label="ANEOS unit tests.")
# Unit tests of the ANEOS equation of state.  This is a silly test, just checking
# that our convenient constructors using the provided inputs build properly.
import unittest
from Spheral1d import *

#===============================================================================
# Unit tests.
#===============================================================================
class TestANEOS(unittest.TestCase):

    #===========================================================================
    # quartz
    #===========================================================================
    def test_quartz(self):
        units = CGS()
        eos = ANEOS(material = "quartz",
                    constants = units)

    #===========================================================================
    # dunite
    #===========================================================================
    def test_dunite(self):
        units = CGS()
        eos = ANEOS(material = "dunite",
                    constants = units)

    #===========================================================================
    # serpentine
    #===========================================================================
    # def test_serpentine(self):
    #     units = CGS()
    #     eos = ANEOS(material = "serpentine",
    #                 constants = units)

#===============================================================================
# Run the suckers.
#===============================================================================
if __name__ == "__main__":
    unittest.main()
