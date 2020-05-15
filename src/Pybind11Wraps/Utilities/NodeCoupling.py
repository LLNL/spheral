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
                 nodeListi = "const unsigned",
                 i = "const unsigned",
                 nodeListj = "const unsigned",
                 j = "const unsigned"):
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
               damageGradient = "const FieldList<%(Dimension)s, Vector>&",
               H = "const FieldList<%(Dimension)s, SymTensor>&"):
        "Constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 nodeListi = "const unsigned",
                 i = "const unsigned",
                 nodeListj = "const unsigned",
                 j = "const unsigned"):
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
               damageGradient = "const FieldList<%(Dimension)s, Vector>&",
               H = "const FieldList<%(Dimension)s, SymTensor>&",
               fragIDs = "const FieldList<%(Dimension)s, int>&"):
        "Constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 nodeListi = "const unsigned",
                 i = "const unsigned",
                 nodeListj = "const unsigned",
                 j = "const unsigned"):
        "Provides a damaged coupling between nodes (nodeListi, i) <-> (nodeListj, j)"
        return "double"
