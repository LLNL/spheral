#-------------------------------------------------------------------------------
# MortonOrderRedistributeNodes
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SpaceFillingCurveRedistributeNodes import *

@PYB11template("Dimension")
class MortonOrderRedistributeNodes(SpaceFillingCurveRedistributeNodes):
    """MortonOrderRedistributeNodes

Attempt to redistribute nodes such that they are laid out in memory
in a Morton ordering.  Note that this involves renumbering the nodes of 
each NodeList, not just redistributing them between processors.

Warren & Salmon (1995), Computer Physics Communications, 87, 266-290."""

    PYB11typedefs = """
    typedef typename KeyTraits::Key Key;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dummy = "const double",
               minNodesPerDomainFraction = ("const double", "0.5"),
               maxNodesPerDomainFraction = ("const double", "1.5"),
               workBalance = ("const bool", "true"),
               localReorderOnly = ("const bool", "false")):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def computeHashedIndices(self,
                             dataBase = "const DataBase<%(Dimension)s>&"):
        "This is the required method for all descendant classes."
        return "FieldList<%(Dimension)s, Key> "
