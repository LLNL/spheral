#-------------------------------------------------------------------------------
# CellFaceFlag
#-------------------------------------------------------------------------------
from PYB11Generator import *

class CellFaceFlag:
    "A struct for encoding cell face information returned by computeVoronoiVolume"

    cellFace = PYB11readwrite(doc="The index of the face in the cell.facets array")
    nodeListj = PYB11readwrite(doc="The NodeList of the opposite node across the face")
    j = PYB11readwrite(doc="The index of the opposite node across the face")

    def pyinit(self):
        "Default constructor"
        return

    def pyinit1(self, x="int", y="int", z="int"):
        "Construct with (cellFace=x, nodeListj=y, j=z)"
        return

    def __eq__(self):
        return

    @PYB11implementation("[](const CellFaceFlag& self) { std::stringstream os; os << self; return os.str(); }")
    def __repr__(self):
        return
