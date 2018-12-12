#-------------------------------------------------------------------------------
# KeyTraits
#-------------------------------------------------------------------------------
from PYB11Generator import *

class KeyTraits:
    """Encapsulate how we think about keys for the space filling curves.
This is specialized assuming we're working with 64 bit uint64_t."""

    numbits =   PYB11readonly(static=True, doc="Total number of bits in a Key")
    numbits1d = PYB11readonly(static=True, doc="Number of bits allowed per dimension in a Key")
    zero =      PYB11readonly(static=True, doc="Null value")
    one =       PYB11readonly(static=True, doc="Unit value")
    two =       PYB11readonly(static=True, doc="Like, the two value")
    maxKey1d =  PYB11readonly(static=True, doc="The maximum value in a dimension")
    maxKey =    PYB11readonly(static=True, doc="The maximum possible key")
