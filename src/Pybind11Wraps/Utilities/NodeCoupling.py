#-------------------------------------------------------------------------------
# NodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11module("SpheralUtilities")
class NodeCoupling:
    "A functor base class encapsulating how we couple pairs of nodes."

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 pair = "const NodePairIdxType&"):
        "Functional method to override for coupling (nodeListi, i) <-> (nodeListj, j)"
        return "double"

#-------------------------------------------------------------------------------
# DamagedNodeCoupling
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class DamagedNodeCoupling(NodeCoupling):
    "A functor class encapsulating how we couple solid nodes in the presence of multiple materials and damage."

    PYB11typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self,
               damage = "const FieldList<%(Dimension)s, SymTensor>&",
               pairs = "NodePairList&"):
        "Constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 pair = "const NodePairIdxType&"):
        "Provides a damaged coupling between nodes (nodeListi, i) <-> (nodeListj, j)"
        return "double"

#-------------------------------------------------------------------------------
# DamagedNodeCouplingWithFrags
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class DamagedNodeCouplingWithFrags(DamagedNodeCoupling):
    """A functor class encapsulating how we couple solid nodes in the presence of
multiple materials and damage.  This version adds logic to decouple based
on fragment ID as well."""

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    def pyinit(self,
               damage = "const FieldList<%(Dimension)s, SymTensor>&",
               fragIDs = "const FieldList<%(Dimension)s, int>&",
               pairs = "NodePairList&"):
        "Constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 pair = "const NodePairIdxType&"):
        "Provides a damaged coupling between nodes (nodeListi, i) <-> (nodeListj, j)"
        return "double"

#-------------------------------------------------------------------------------
# ThreePointDamagedNodeCoupling
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class ThreePointDamagedNodeCoupling(NodeCoupling):
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
               position = "const FieldList<%(Dimension)s, Vector>&",
               H = "const FieldList<%(Dimension)s, SymTensor>&",
               damage = "const FieldList<%(Dimension)s, SymTensor>&",
               W = "const TableKernel<%(Dimension)s>&",
               connectivity = "const ConnectivityMap<%(Dimension)s>&",
               useIntersectConnectivity = "const bool",
               pairs = "NodePairList&"):
        "Constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 pair = "const NodePairIdxType&"):
        "Provides a damaged coupling between nodes (nodeListi, i) <-> (nodeListj, j)"
        return "double"
