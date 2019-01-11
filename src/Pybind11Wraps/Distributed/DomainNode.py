#-------------------------------------------------------------------------------
# DomainNode
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class DomainNode:

    #...........................................................................
    # Attributes
    localNodeID = PYB11readwrite()
    uniqueLocalNodeID = PYB11readwrite()
    globalNodeID = PYB11readwrite()
    nodeListID = PYB11readwrite()
    domainID = PYB11readwrite()
    work = PYB11readwrite()
    position = PYB11readwrite()

    packSize = PYB11property("size_t", static=True)
    pack = PYB11property("std::vector<double>")

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    
