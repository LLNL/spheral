"""
polytope module.

Rewrap the polytope interface for Python with pybind11.  This should
be superseeded by natively redoing polytope's Python interface with 
PYB11 at some point.
"""

from PYB11Generator import *

PYB11includes = ['"polytope/Tessellation.hh"',
                 '<sstream>']

PYB11namespaces = ["polytope"]

from Tessellation import *

# Instantiate the types
Tessellation2d = PYB11TemplateClass(Tessellation, template_parameters=("2", "double"))
Tessellation3d = PYB11TemplateClass(Tessellation, template_parameters=("3", "double"))
