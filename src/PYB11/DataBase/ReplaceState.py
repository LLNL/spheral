from PYB11Generator import *
from FieldUpdatePolicy import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class ReplaceState(FieldUpdatePolicy):

    PYB11typedefs = """
    using KeyType = typename ReplaceState<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"
        return

    @PYB11implementation("[](py::list depends) { auto result = make_policy<ReplaceState<%(Dimension)s, %(ValueType)s>>(); for (auto x: depends) result->addDependency(x.cast<std::string>()); return result; }")
    def pyinit(self,
               depends = ("py::list", "py::list()")):
        "Build with some set of dependencies"
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
        "Prefix for key for new values"
        return "const std::string"
