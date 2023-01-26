#-------------------------------------------------------------------------------
# ThreePointDamagedNodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *
import NodeCoupling

@PYB11template("Dimension")
class ThreePointDamagedNodeCoupling(NodeCoupling.NodeCoupling):
    """A functor class encapsulating how we couple solid nodes in the presence of
multiple materials and damage.

This for uses the "three point" formalism, which allows damaged points to
cut communication between pairs that talk across them."""

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self,
               state = "const State<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               pairs = "NodePairList&"):
        "Constructor"
