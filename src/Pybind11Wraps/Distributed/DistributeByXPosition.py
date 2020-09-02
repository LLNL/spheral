#-------------------------------------------------------------------------------
# DistributeByXPosition
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RedistributeNodes import *

@PYB11template("Dimension")
class DistributeByXPosition(RedistributeNodes):
    """DistributeByXPosition -- Redistribute nodes by sorting their positions
in x coordinate.  Really only useful in 1-D, as a test."""

    PYB11typedefs = """
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def redistributeNodes(self,
                          dataBase = "DataBase<%(Dimension)s>&",
                          boundaries = ("std::vector<Boundary<%(Dimension)s>*>", "std::vector<Boundary<%(Dimension)s>*>()")):
        """Given a Spheral++ data base of NodeLists, repartition it among the processors.
This is the method required of all descendent classes."""
        return "void"
