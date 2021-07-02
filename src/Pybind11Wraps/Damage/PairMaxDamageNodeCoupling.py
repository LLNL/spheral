#-------------------------------------------------------------------------------
# PairMaxDamageNodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *
import NodeCoupling

@PYB11template("Dimension")
class PairMaxDamageNodeCoupling(NodeCoupling.NodeCoupling):
    """A functor class encapsulating how we couple solid nodes in the presence of
multiple materials and damage.

This form simply directly damages each pair based on their mutual damage."""

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self,
               state = "const State<%(Dimension)s>&",
               pairs = "NodePairList&"):
        "Constructor"

