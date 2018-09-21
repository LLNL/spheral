from PYB11Generator import *
from NodeList import NodeList

#-------------------------------------------------------------------------------
# FluidNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class FluidNodeList(NodeList):
    "Spheral FluidNodeList base class in %(Dimension)s, i.e.,  the NodeList for fluid hydrodynamics."

    def pyinit(self,
               name = "std::string",
               eos = "EquationOfState<%(Dimension)s>&",
               numInternal = ("unsigned", "0"),
               numGhost = ("unsigned", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("unsigned", "500"),
               rhoMin = ("double", "1e-10"),
               rhoMax = ("double", "1e100")):
        "Constructor for a FluidNodeList class."
        return

    @PYB11const
    def massDensity(self):
        "The mass density field"
        return "const %(ScalarField)s&"

    @PYB11const
    def specificThermalEnergy(self):
        "The specific thermal energy field"
        return "const %(ScalarField)s&"

    @PYB11virtual
    @PYB11const
    def pressure(self, result="%(ScalarField)s&"):
        "Compute the current pressure, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11virtual
    @PYB11const
    def temperature(self, result="%(ScalarField)s&"):
        "Compute the current temperature, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11virtual
    @PYB11const
    def soundSpeed(self, result="%(ScalarField)s&"):
        "Compute the current sound speed, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11virtual
    @PYB11const
    def volume(self, result="%(ScalarField)s&"):
        "Compute the current volume, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11virtual
    @PYB11const
    def linearMomentum(self, result="%(VectorField)s&"):
        "Compute the current linear momentum, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11virtual
    @PYB11const
    def totalEnergy(self, result="%(ScalarField)s&"):
        "Compute the current total energy, storing the result in the argument %(ScalarField)s."
        return "void"

    @PYB11const
    def equationOfState(self):
        "Return the equation of state object this FluidNodeList is associated with."
        return "const EquationOfState<%(Dimension)s>&"

    def setequatinoOfState(self, equationOfState="const EquationOfState<%(Dimension)s>&"):
        "Set the equation of state for this FluidNodeList."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"

    # Comparison
    def __eq__(self):
        "Equivalence test with another FluidNodeList"

    def __ne__(self):
        "Inequivalence test with another FluidNodeList"

    # Methods used for properties
    @PYB11ignore
    @PYB11cppname("rhoMin")
    @PYB11const
    def getrhoMin(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("rhoMin")
    @PYB11const
    def setrhoMin(self, val="double"):
        return "void"

    @PYB11ignore
    @PYB11cppname("rhoMax")
    @PYB11const
    def getrhoMax(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("rhoMax")
    @PYB11const
    def setrhoMax(self, val="double"):
        return "void"

    # Properties
    rhoMin = property(getrhoMin, setrhoMin, doc="The minimum allowed mass density.")
    rhoMax = property(getrhoMax, setrhoMax, doc="The maximum allowed mass density.")
