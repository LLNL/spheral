#-------------------------------------------------------------------------------
# uniform_random_01
#-------------------------------------------------------------------------------
from PYB11Generator import *

class uniform_random_01:
    "Encapsulate a random number generator to generate numbers in [0,1)."

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit1(self, seed="const size_t"):
        "Construct with a seed value"

    #...........................................................................
    # Methods
    def __call__(self):
        "Return a random value in the range [0,1)"
        return "double"

    def seed(self, val="const size_t"):
        "Set the seed for our generator"
        return "void"

    def advance(self, n="const size_t"):
        "Advance the state of the generator as though called n times"
        return "void"

    def serialize(self, buffer="std::vector<char>"):
        "Serialize to a buffer of bytes"
        return "void"

    @PYB11implementation("""
    [](uniform_random_01& self, std::vector<char>& buffer, size_t pos) {
      auto itr = buffer.begin() + pos;
      self.deserialize(itr, buffer.end());
      return std::distance(buffer.begin(), itr);
    }""")
    def deserialize(self,
                    buffer = "std::vector<char>",
                    pos = "size_t"):
        "Deserialize from a buffer of bytes, starting at the given position.  Returns the new position in the buffer."
        return "size_t"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    def __ne__(self):
        return
