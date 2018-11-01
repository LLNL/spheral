"""
Spheral NodeList module.

Provides the fundamental NodeList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
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
            '"Material/EquationOfState.hh"',
            '"SolidMaterial/StrengthModel.hh"',
            '"Kernel/TableKernel.hh"',
            '"Neighbor/ConnectivityMap.hh"',
            '"Mesh/Mesh.hh"',
            '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# preamble
#-------------------------------------------------------------------------------
for ndim in dims:
    preamble += "typedef std::pair<NodeList<Dim<%(ndim)i>>*, std::string> pair_NodeList%(ndim)idptr_string;\n" % {"ndim": ndim}

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
from FluidNodeList import FluidNodeList
from SolidNodeList import SolidNodeList
from SmoothingScaleBase import SmoothingScaleBase
from FixedSmoothingScale import FixedSmoothingScale
from SPHSmoothingScale import SPHSmoothingScale
from ASPHSmoothingScale import ASPHSmoothingScale

for ndim in dims:
    exec('''
NodeListRegistrar%(ndim)id = PYB11TemplateClass(NodeListRegistrar, template_parameters="Dim<%(ndim)i>")

NodeList%(ndim)id = PYB11TemplateClass(NodeList, template_parameters="Dim<%(ndim)i>")
FluidNodeList%(ndim)id = PYB11TemplateClass(FluidNodeList, template_parameters="Dim<%(ndim)i>")
SolidNodeList%(ndim)id = PYB11TemplateClass(SolidNodeList, template_parameters="Dim<%(ndim)i>")

SmoothingScaleBase%(ndim)id = PYB11TemplateClass(SmoothingScaleBase, template_parameters="Dim<%(ndim)i>")
FixedSmoothingScale%(ndim)id = PYB11TemplateClass(FixedSmoothingScale, template_parameters="Dim<%(ndim)i>")
SPHSmoothingScale%(ndim)id = PYB11TemplateClass(SPHSmoothingScale, template_parameters="Dim<%(ndim)i>")
ASPHSmoothingScale%(ndim)id = PYB11TemplateClass(ASPHSmoothingScale, template_parameters="Dim<%(ndim)i>")

vector_of_NodeList%(ndim)id = PYB11_bind_vector("NodeList<Dim<%(ndim)i>>*", opaque=True, local=False)
vector_of_FluidNodeList%(ndim)id = PYB11_bind_vector("FluidNodeList<Dim<%(ndim)i>>*", opaque=True, local=False)
vector_of_SolidNodeList%(ndim)id = PYB11_bind_vector("SolidNodeList<Dim<%(ndim)i>>*", opaque=True, local=False)

vector_of_pair_NodeList%(ndim)id_string = PYB11_bind_vector("pair_NodeList%(ndim)idptr_string", opaque=True, local=False)
''' % {"ndim" : ndim})
