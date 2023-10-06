from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class FieldUpdatePolicy(UpdatePolicyBase):
    "FieldUpdatePolicy -- Base/interface class for the policies on how Field state variables are to be updated."

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
    @PYB11const
    def clonePerField(self):
        "Returns whether this policy should be cloned for each Field in a FieldList or not"
        return "bool"
