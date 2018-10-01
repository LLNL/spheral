#-------------------------------------------------------------------------------
# AxisBoundaryRZ
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *

@PYB11template()
class AxisBoundaryRZ(Boundary):

    typedefs = """
    typedef Dim<2> DIM;
    typedef DIM::Scalar Scalar;
    typedef DIM::Vector Vector;
    typedef DIM::Tensor Tensor;
    typedef DIM::SymTensor SymTensor;
    typedef DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               etamin = "const double"):
        "Construct with the DataBase"

    #...........................................................................
    # Methods
    @PYB11virtual
    def setViolationNodes(self, nodeList="NodeList<DIM>&"):
        return "void"

    @PYB11virtual
    def updateViolationNodes(self, nodeList="NodeList<DIM>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "The label for writing in restart files"
        return "std::string"

    #...........................................................................
    # Properties
    etamin = PYB11property("const double", "etamin", "etamin", doc="The fuzz value for approaching the axis")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
#PYB11inject(BoundaryAbstractMethods, AxisBoundaryRZ, virtual=True, pure_virtual=False)
