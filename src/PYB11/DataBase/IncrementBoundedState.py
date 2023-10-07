from PYB11Generator import *
from FieldUpdatePolicy import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class IncrementBoundedState(FieldUpdatePolicy):

    PYB11typedefs = """
    using KeyType = typename IncrementBoundedState<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with no dependencies"
        return

    def pyinit1(self,
                depend0 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with one dependencies"
        return

    def pyinit2(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with two dependencies"
        return

    def pyinit3(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with three dependencies"
        return

    def pyinit4(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with four dependencies"
        return

    def pyinit5(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with five dependencies"
        return

    def pyinit6(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                depend5 = "const std::string&",
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
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
        "Update a Field assoicated with the given key"
        return "void"

    @PYB11static
    def prefix(self):
        "Prefix for key of derivatives"
        return "const std::string"

    #...........................................................................
    # Properties
    minValue = PYB11property(doc="Minimum bound for Field values")
    maxValue = PYB11property(doc="Maximum bound for Field values")
