
from PYB11Generator import *
#from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes = ['"Geometry/Dimension.hh"',
                 '"axom/sidre.hpp"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

# from Group import *
from View import *

#-------------------------------------------------------------------------------
# Instantiate types and add dimension dependent functions.
#-------------------------------------------------------------------------------
