"""
Spheral Neighbor module.

Provides the base class and implementations for neighbor finding in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Boundary/Boundary.hh"',
                  '"Neighbor/GridCellIndex.hh"',
                  '"Neighbor/GridCellPlane.hh"',
                  '"Neighbor/Neighbor.hh"',
                  '"Neighbor/NestedGridNeighbor.hh"',
                  '"Neighbor/TreeNeighbor.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"Neighbor/NodePairList.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
NeighborSearchType = PYB11enum(("Gather", "Scatter", "GatherScatter"), export_values=True,
                               doc="Choice for how nodes should look for neighbors")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from GridCellIndex import *
from GridCellPlane import *
from Neighbor import *
from NestedGridNeighbor import *
from TreeNeighbor import *
from ConnectivityMap import *
from NodePairList import *

for ndim in dims:
    exec('''
GridCellIndex%(ndim)id = PYB11TemplateClass(GridCellIndex, template_parameters="Dim<%(ndim)i>")
GridCellPlane%(ndim)id = PYB11TemplateClass(GridCellPlane, template_parameters="Dim<%(ndim)i>")
Neighbor%(ndim)id = PYB11TemplateClass(Neighbor, template_parameters="Dim<%(ndim)i>")
NestedGridNeighbor%(ndim)id = PYB11TemplateClass(NestedGridNeighbor, template_parameters="Dim<%(ndim)i>")
TreeNeighbor%(ndim)id = PYB11TemplateClass(TreeNeighbor, template_parameters="Dim<%(ndim)i>")
ConnectivityMap%(ndim)id = PYB11TemplateClass(ConnectivityMap, template_parameters="Dim<%(ndim)i>")

vector_of_GridCellIndex%(ndim)id = PYB11_bind_vector("GridCellIndex<Dim<%(ndim)i>>", opaque=True, local=False)
vector_of_vector_of_GridCellIndex%(ndim)id = PYB11_bind_vector("std::vector<GridCellIndex<Dim<%(ndim)i>>>", opaque=True, local=False)
''' % {"ndim"      : ndim})
