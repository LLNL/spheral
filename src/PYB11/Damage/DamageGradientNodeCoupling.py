#-------------------------------------------------------------------------------
# DamageGradientNodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *
import NodeCoupling

@PYB11template("Dimension")
class DamageGradientNodeCoupling(NodeCoupling.NodeCoupling):
    """A functor class encapsulating how we couple solid nodes in the presence of
multiple materials and damage.

This one attempts to mock up the shielding effect of ThreePointDamagedNodeCoupling
by using local damage gradient to estimate when nodes are separated by
regions of greater damage (or fractures)."""

