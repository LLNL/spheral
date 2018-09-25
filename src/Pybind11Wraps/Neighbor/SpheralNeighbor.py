"""
Spheral Neighbor module.

Provides the base class and implementations for neighbor finding in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Geometry/GeomPlane.hh"',
            '"Neighbor/GridCellIndex.hh"',
            '"Neighbor/GridCellPlane.hh"',
            '"Neighbor/Neighbor.hh"',
            '"Neighbor/NestedGridNeighbor.hh"',
            '"Neighbor/TreeNeighbor.hh"',
            '"Neighbor/ConnectivityMap.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Define a preamble function to expose the protected methods of Neighbor.
#-------------------------------------------------------------------------------
preamble = """
namespace Spheral {
  template<typename Dimension>
  class NeighborPublicist: public Neighbor<Dimension> {
  public:
    using Neighbor<Dimension>::accessNodeExtentField;
  };
}
"""

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
NeighborSearchType = PYB11enum(("None", "Gather", "Scatter", "GatherScatter"), export_values=True,
                               doc="Choice for how nodes should look for neighbors")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from GridCellIndex import *
from Neighbor import *

for ndim in dims:
    exec('''
GridCellIndex%(ndim)id = PYB11TemplateClass(GridCellIndex, template_parameters="Dim<%(ndim)i>")
Neighbor%(ndim)id = PYB11TemplateClass(Neighbor, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})
