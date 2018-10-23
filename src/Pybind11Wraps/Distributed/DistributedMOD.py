"""
Spheral Boundary module.

Provides the Boundary base class and many boundary implementations.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

from DistributedBoundary import *
from NestedGridDistributedBoundary import *
from TreeDistributedBoundary import *
from BoundingVolumeDistributedBoundary import *
from DomainNode import *
from RedistributeNodes import *
from SpaceFillingCurveRedistributeNodes import *
from MortonOrderRedistributeNodes import *
from PeanoHilbertOrderRedistributeNodes import *
from SortAndDivideRedistributeNodes import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Boundary/Boundary.hh"',
            '"Neighbor/GridCellIndex.hh"',
            '"Distributed/DistributedBoundary.hh"',
            '"Distributed/NestedGridDistributedBoundary.hh"',
            '"Distributed/BoundingVolumeDistributedBoundary.hh"',
            '"Distributed/TreeDistributedBoundary.hh"',
            '"Distributed/DomainNode.hh"',
            '"Distributed/RedistributeNodes.hh"',
            '"Distributed/DistributeByXPosition.hh"',
            '"Distributed/SpaceFillingCurveRedistributeNodes.hh"',
            '"Distributed/MortonOrderRedistributeNodes.hh"',
            '"Distributed/PeanoHilbertOrderRedistributeNodes.hh"',
            '"Distributed/SortAndDivideRedistributeNodes1d.hh"',
            '"Distributed/SortAndDivideRedistributeNodes2d.hh"',
            '"Distributed/SortAndDivideRedistributeNodes3d.hh"',
            '"Distributed/VoronoiRedistributeNodes.hh"',
            '"Utilities/KeyTraits.hh"',
            '"Field/Field.hh"',
            '"Field/FieldList.hh"',
            '"FileIO/FileIO.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
DistributedBoundary%(ndim)id = PYB11TemplateClass(DistributedBoundary, template_parameters="%(Dimension)s")
NestedGridDistributedBoundary%(ndim)id = PYB11TemplateClass(NestedGridDistributedBoundary, template_parameters="%(Dimension)s")
TreeDistributedBoundary%(ndim)id = PYB11TemplateClass(TreeDistributedBoundary, template_parameters="%(Dimension)s")
BoundingVolumeDistributedBoundary%(ndim)id = PYB11TemplateClass(BoundingVolumeDistributedBoundary, template_parameters="%(Dimension)s")
DomainNode%(ndim)id = PYB11TemplateClass(DomainNode, template_parameters="%(Dimension)s")
RedistributeNodes%(ndim)id = PYB11TemplateClass(RedistributeNodes, template_parameters="%(Dimension)s")
SpaceFillngCurveRedistributeNodes%(ndim)id = PYB11TemplateClass(SpaceFillingCurveRedistributeNodes, template_parameters="%(Dimension)s")
MortonOrderRedistributeNodes%(ndim)id = PYB11TemplateClass(MortonOrderRedistributeNodes, template_parameters="%(Dimension)s")
PeanoHilbertOrderRedistributeNodes%(ndim)id = PYB11TemplateClass(PeanoHilbertOrderRedistributeNodes, template_parameters="%(Dimension)s")
SortAndDivideRedistributeNodes%(ndim)id = PYB11TemplateClass(SortAndDivideRedistributeNodes, template_parameters="%(Dimension)s")

#vector_of_Boundary%(ndim)id = PYB11_bind_vector("Boundary<%(Dimension)s>*", opaque=True, local=False)
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})

if 1 in dims:
    from SortAndDivideRedistributeNodes1d import *

if 2 in dims:
    from SortAndDivideRedistributeNodes2d import *

if 3 in dims:
    from SortAndDivideRedistributeNodes3d import *
