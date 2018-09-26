"""
Spheral FileIO module.

Provides the interfaces for communicating with file systems.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"FileIO/FileIO.hh"',
            '"FileIO/FlatFileIO.hh"',
            '"FileIO/SiloFileIO.hh"',
            '"FileIO/PyFileIO.hh"',
            '"FileIO/vectorstringUtilities.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
AccessType = PYB11enum(("Undefined", 
                        "Create", 
                        "Read",
                        "Write",
                        "ReadWrite"), export_values=True,
                       doc="How are we opening/accessing a file")

FlatFileFormat = PYB11enum(("ascii", "binary"), export_values=True,
                           doc="Format of ascii file")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from FileIO import *

# for ndim in dims:
#     exec('''
# FileIO%(ndim)id = PYB11TemplateClass(FileIO, template_parameters="Dim<%(ndim)i>")
# ''' % {"ndim" : ndim})
