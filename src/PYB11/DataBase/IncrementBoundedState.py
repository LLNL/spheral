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
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())"),
                wildCardDerivs = ("const bool", "false")):
        "Build with no dependencies"
        return

    @PYB11implementation("[](py::list depends, const %(ValueType)s minValue, const %(ValueType)s maxValue, const bool wildCardDerivs) { auto result = make_policy<IncrementBoundedState<%(Dimension)s, %(ValueType)s>>(minValue, maxValue, wildCardDerivs); for (auto x: depends) result->addDependency(x.cast<std::string>()); return result; }")
    def pyinit1(self,
                depends = ("py::list", "py::list()"),
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with dependencies"
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
    wildCardDerivs = PYB11property("bool", getter="wildCardDerivs", setter="wildCardDerivs")
