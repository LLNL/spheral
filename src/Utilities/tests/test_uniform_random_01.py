#ATS:test(SELF, label="Unit tests of uniform_random_01 random number generator")
from math import *
from Spheral import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
rangen = random.Random()

ntests = 1000

class TestRandom01(unittest.TestCase):

    #===========================================================================
    # Various ways of constructing
    #===========================================================================
    def testConstructors(self):
        seed = rangen.randint(1, 2**64)
        gen1 = uniform_random_01()
        gen1.seed(seed)
        gen2 = uniform_random_01(seed)
        gen3 = uniform_random_01(gen1)
        assert gen1 == gen2
        assert gen1 == gen3
        for i in xrange(ntests):
            assert gen1() == gen2() == gen3()

    #===========================================================================
    # Comparisons
    #===========================================================================
    def testComparisons(self):
        seed = rangen.randint(1, 2**64)
        gen1 = uniform_random_01(seed)
        gen2 = uniform_random_01()
        assert gen1 != gen2
        gen2.seed(seed)
        assert gen1 == gen2

    #===========================================================================
    # advance
    #===========================================================================
    def testAdvance(self):
        seed = rangen.randint(1, 2**64)
        gen1 = uniform_random_01(seed)
        throwaway = [gen1() for i in xrange(ntests)]
        vals1 = [gen1() for i in xrange(ntests)]
        gen2 = uniform_random_01(seed)
        gen2.advance(ntests)
        vals2 = [gen2() for i in xrange(ntests)]
        assert vals1 == vals2

    #===========================================================================
    # Serialization
    #===========================================================================
    def testSerialize(self):
        seed = rangen.randint(1, 2**64)
        gen1 = uniform_random_01(seed)
        throwaway = [gen1() for i in xrange(ntests)]
        buf = vector_of_char()
        gen1.serialize(buf)
        gen2 = uniform_random_01()
        i = gen2.deserialize(buf, 0)
        assert i == len(buf)
        vals1 = [gen1() for i in xrange(ntests)]
        vals2 = [gen2() for i in xrange(ntests)]
        assert vals1 == vals2

if __name__ == "__main__":
    unittest.main()

