from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class NodeList:
    "Spheral NodeList base class in %(Dimension)s"

    def pyinit(self,
               name = "std::string",
               numInternal = ("unsigned", "0"),
               numGhost = ("unsigned", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("unsigned", "500")):
        "Constructor for NodeList base class."
        return

    def nodeType(self,
                 ID = "int"):
        "Return the classification of the given node"
        return "NodeType"

    @PYB11const
    def mass(self):
        "The mass field"
        return "const Field<%(Dimension)s, %(Scalar)s>&"

    @PYB11const
    def positions(self):
        "The position field"
        return "const Field<%(Dimension)s, %(Vector)s>&"

    @PYB11const
    def velocity(self):
        "The velocity field"
        return "const Field<%(Dimension)s, %(Vector)s>&"

    @PYB11const
    def Hfield(self):
        "The H tensor field"
        return "const Field<%(Dimension)s, %(SymTensor)s>&"

    @PYB11const
    def work(self):
        "The CPU work field"
        return "Field<%(Dimension)s, %(Scalar)s>&"

    @PYB11pycppname("mass")
    def setmass(self, newValue="const Field<%(Dimension)s, %(Scalar)s>&"):
        "Set the mass field"
        return "void"

    @PYB11pycppname("positions")
    def setpositions(self, newValue="const Field<%(Dimension)s, %(Vector)s>&"):
        "Set the position field"
        return "void"

    @PYB11pycppname("velocity")
    def setvelocity(self, newValue="const Field<%(Dimension)s, %(Vector)s>&"):
        "Set the velocity field"
        return "void"

    @PYB11pycppname("Hfield")
    def setHfield(self, newValue="const Field<%(Dimension)s, %(SymTensor)s>&"):
        "Set the H tensor field"
        return "void"

    @PYB11pycppname("work")
    def setwork(self, newValue="const Field<%(Dimension)s, %(Scalar)s>&"):
        "Set the CPU work field"
        return "void"
