#-------------------------------------------------------------------------------
# CRKSPHVoidBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class CRKSPHVoidBoundary(Boundary):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef GeomPlane<DIM> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               surfacePoint = "const FieldList<DIM, int>&",
               etaVoidPoints = "const FieldList<DIM, std::vector<Vector>>&"):
        "Constructor"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, CRKSPHVoidBoundary, virtual=True, pure_virtual=False)
