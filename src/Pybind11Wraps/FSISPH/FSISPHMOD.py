"""
Spheral SPH module.

Provides implementations of SPH, PSPH, and ASPH 
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FSISPHHydroBase import *
from SolidFSISPHHydroBase import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"FSISPH/SolidFSISPHHydroBase.hh"',
                  '"FSISPH/SolidFSISPHHydroBaseRZ.hh"',
                  '"FSISPH/FSISPHHydroBase.hh"',
                  '"FSISPH/FSISPHHydroBaseRZ.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
SolidFSISPHHydroBase%(ndim)id = PYB11TemplateClass(SolidFSISPHHydroBase, template_parameters="%(Dimension)s")
FSISPHHydroBase%(ndim)id = PYB11TemplateClass(FSISPHHydroBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

if 2 in dims:
    from FSISPHHydroBaseRZ import *
    from SolidFSISPHHydroBaseRZ import *
