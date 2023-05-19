"""
Spheral ANEOS module.

Provides the ANEOS equation of state
"""

from PYB11Generator import *
#from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from ANEOS import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes = ['"Geometry/Dimension.hh"',
                 '"Material/PhysicalConstants.hh"',
                 '"Material/EquationOfState.hh"',
                 '"Field/Field.hh"',
                 '"SolidMaterial/ANEOS.hh"',
                 '"ANEOSWrappers.hh"',
                 '<vector>',
                 '<string>',
                 '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
def initializeANEOS(in_filename = "std::string",
                    out_filename = "std::string",
                    izetl = "std::vector<int>"):
    """Initialize ANEOS with some rather arcane input.

in_filename  : The name of the ANEOS input file, initializes the Fortran ANEOS library
out_filename : An optional file to write any output from the ANEOS intialization call
izetl        : An array of the material numbers ("EOS#" in the ANEOS input file)
               Note, these numbers must be the negative of the "EOS#" in the input
"""
    return "void"

#-------------------------------------------------------------------------------
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
ANEOS%(ndim)id = PYB11TemplateClass(ANEOS, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
