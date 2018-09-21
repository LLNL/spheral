"""
Spheral NodeList module.

Provides the fundamental NodeList classes.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"NodeList/NodeListRegistrar.hh"',
            '"NodeList/NodeList.hh"',
            '"NodeList/FluidNodeList.hh"',
            '"NodeList/SolidNodeList.hh"',
            '"NodeList/SmoothingScaleBase.hh"',
            '"NodeList/FixedSmoothingScale.hh"',
            '"NodeList/SPHSmoothingScale.hh"',
            '"NodeList/ASPHSmoothingScale.hh"',
            '"NodeList/generateVoidNodes.hh"',
            '"NodeList/nthNodalMoment.hh"',
            '"Kernel/TableKernel.hh"',
            '"Mesh/Mesh.hh"',
            '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
NodeType = PYB11enum(("InternalNode", "GhostNode"), export_values=True,
                     doc="The classifications of Spheral nodes.")

#-------------------------------------------------------------------------------
# NodeListRegistrar
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11singleton
class NodeListRegistrar:

    @PYB11const
    def valid(self):
        return "bool"

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instance")
    @PYB11ignore
    def getinstance(self):
        return "NodeListRegistrar<%(Dimension)s>&"
    instance = property(getinstance, doc="The static NodeListRegistrar<%(Dimension)s> instance.")

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from NodeList import NodeList

for ndim in dims:
    exec('''
NodeListRegistrar%(ndim)id = PYB11TemplateClass(NodeListRegistrar, template_parameters=dimDictionary(%(ndim)i))
NodeList%(ndim)id = PYB11TemplateClass(NodeList, template_parameters = dimDictionary(%(ndim)i))
''' % {"ndim" : ndim})

