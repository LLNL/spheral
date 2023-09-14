from PYB11Generator import *
from FieldListUpdatePolicyBase import *

@PYB11template("Dimension", "ValueType")
class IncrementFieldList(FieldListUpdatePolicyBase):

    PYB11typedefs = """
    using KeyType = typename IncrementFieldList<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                wildCardDerivs = ("const bool", "false")):
        "Build with no dependencies"
        return

    def pyinit1(self,
                depend0 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with one dependencies"
        return

    def pyinit2(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with two dependencies"
        return

    def pyinit3(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with three dependencies"
        return

    def pyinit4(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with four dependencies"
        return

    def pyinit5(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with five dependencies"
        return

    def pyinit6(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                depend5 = "const std::string&",
                wildCardDerivs = ("const bool", "false")):
        "Build with six dependencies"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    def update(self,
               key = "const KeyType&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double"):
        "Update a FieldList assoicated with the given key"
        return "void"

    @PYB11static
    def prefix(self):
        "Prefix for key of derivatives"
        return "const std::string"
