from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("DataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class UpdatePolicyBase:

    PYB11typedefs = """
    using KeyType = typename UpdatePolicyBase<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Build with no dependencies"
        return

    def pyinit1(self,
                depend0 = "const std::string&"):
        "Build with one dependencies"
        return

    def pyinit2(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&"):
        "Build with two dependencies"
        return

    def pyinit3(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&"):
        "Build with three dependencies"
        return

    def pyinit4(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&"):
        "Build with four dependencies"
        return

    def pyinit5(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&"):
        "Build with five dependencies"
        return

    def pyinit6(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                depend5 = "const std::string&"):
        "Build with six dependencies"
        return

    def pyinit7(self,
                depend0 = "const std::string&",
                depend1 = "const std::string&",
                depend2 = "const std::string&",
                depend3 = "const std::string&",
                depend4 = "const std::string&",
                depend5 = "const std::string&",
                depend6 = "const std::string&"):
        "Build with seven dependencies"
        return

    #...........................................................................
    # Virtual methods
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

    @PYB11virtual
    def updateAsIncrement(self,
                          key = "const KeyType&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&",
                          multiplier = "const double",
                          t = "const double",
                          dt = "const double"):
        """An alternate method to be called when you want to specify that the derivative information
should be assumed to not necessarily be properly time-centered, and therefore you should 
only use time advancement ideas, no "replace" or more sophisticated approaches.
Default to just calling the generic method."""
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def independent(self):
        return "bool"

    @PYB11const
    def dependent(self):
        return "bool"

    @PYB11const
    def dependencies(self):
        "Return the set of field names that this state depends upon (if any)."
        return "const std::vector<std::string>&"

    def addDependency(self,
                      depend = "const std::string&"):
        "Allow the addition of new dependencies."
        return "void"

    @PYB11static
    def wildcard(self):
        "The wildcard string for comparing dependencies"
        return "std::string"

    #...........................................................................
    # Comparators
    @PYB11pure_virtual
    @PYB11const
    @PYB11cppname("operator==")
    @PYB11pyname("equals")      # Prevent converting this to a normal python comparison since it's a pure virtual function
    def __eq__(self,
               rhs = "const UpdatePolicyBase<%(Dimension)s>&"):
        return "bool"
