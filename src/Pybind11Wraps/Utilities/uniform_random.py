#-------------------------------------------------------------------------------
# uniform_random
#-------------------------------------------------------------------------------
from PYB11Generator import *

class uniform_random:
    "Encapsulate a random number generator to generate numbers in [0,1)."

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                seed = "const size_t",
                minVal = "const double",
                maxVal = "const double"):
        "Construct with a seed value and range"

    def pyinit2(self, rhs="const uniform_random&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    def __call__(self):
        "Return a random value in the range [0,1)"
        return "double"

    def advance(self, n="const size_t"):
        "Advance the state of the generator as though called n times"
        return "void"

    def range(self,
              a = "const double",
              b = "const double"):
        "Set the range of generated values to be [a,b)"
        return "void"

    def serialize(self, buffer="std::vector<char>&"):
        "Serialize to a buffer of bytes"
        return "void"

    @PYB11implementation("""
    [](uniform_random& self, const std::vector<char>& buffer, size_t pos) -> size_t {
      auto itr = buffer.begin() + pos;
      self.deserialize(itr, buffer.end());
      return std::distance(buffer.begin(), itr);
    }""")
    def deserialize(self,
                    buffer = "const std::vector<char>&",
                    pos = "size_t"):
        "Deserialize from a buffer of bytes, starting at the given position.  Returns the new position in the buffer."
        return "size_t"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    def __ne__(self):
        return

    #...........................................................................
    # Properties
    seed = PYB11property("size_t", "seed", "seed", doc="Set/get the random number seed")
    numCalls = PYB11property("size_t", "numCalls",
                             doc = "Number of times a random number has been generated")
    min = PYB11property("double", "min",
                        doc = "The minimum number in the generated range")
    max = PYB11property("double", "max",
                        doc = "The maximum number in the generated range")
    
