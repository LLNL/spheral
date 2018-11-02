"""
Spheral DataBase module.

Provides higher level interface for State and DataBase objects.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# More preamble
#-------------------------------------------------------------------------------
preamble += """
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<3>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<3>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<3>>*>)
"""

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"DataBase/DataBase.hh"',
            '"DataBase/StateBase.hh"',
            '"DataBase/State.hh"',
            '"DataBase/StateDerivatives.hh"',
            '"Field/Field.hh"',
            '"Neighbor/ConnectivityMap.hh"',
            '"Physics/Physics.hh"',
            '"Utilities/DataTypeTraits.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from StateBase import *
from State import *
from StateDerivatives import *
from DataBase import *

for ndim in dims:
    exec('''
StateBase%(ndim)id = PYB11TemplateClass(StateBase, template_parameters="Dim<%(ndim)i>")
State%(ndim)id = PYB11TemplateClass(State, template_parameters="Dim<%(ndim)i>")
StateDerivatives%(ndim)id = PYB11TemplateClass(StateDerivatives, template_parameters="Dim<%(ndim)i>")
DataBase%(ndim)id = PYB11TemplateClass(DataBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})
