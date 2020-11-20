
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

from Group import *
from View import *

#-------------------------------------------------------------------------------
# Instantiate types and add dimension dependent functions.
#-------------------------------------------------------------------------------

# @PYB11namespace("axom::sidre::View")
# @PYB11implementation("""[](axom::sidre::View inputView, int n) { 
#                                             int* viewData = inputView.getData();
#                                             py::list result;
#                                             for (int i; i < n; i++) 
#                                                 result.append(viewData[i]);
#                                             return result;
#                                            }""")
# def getDataA(inputView = "axom::sidre::View", n = "int"):
#     "Return data held by view and cast it to any compatible type allowed by Conduit (return type depends on type caller assigns it to)."
#     return "axom::sidre::Node::Value"