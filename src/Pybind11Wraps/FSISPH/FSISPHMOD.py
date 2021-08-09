"""
Spheral FSISPH module.

"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()


from SolidFSISPHHydroBase import *
from SlideSurface import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"FSISPH/SolidFSISPHHydroBase.hh"',
                  '"FSISPH/SlideSurface.hh"',
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
SlideSurface%(ndim)id = PYB11TemplateClass(SlideSurface, template_parameters="%(Dimension)s")
SolidFSISPHHydroBase%(ndim)id = PYB11TemplateClass(SolidFSISPHHydroBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

