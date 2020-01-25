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
PYB11includes += ['"NodeList/NodeListRegistrar.hh"',
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
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# preamble
#-------------------------------------------------------------------------------
for ndim in dims:
    PYB11preamble += "typedef std::pair<NodeList<Dim<%(ndim)i>>*, std::string> pair_NodeList%(ndim)idptr_string;\n" % {"ndim": ndim}

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
NodeType = PYB11enum(("InternalNode", "GhostNode"), export_values=True,
                     doc="The classifications of Spheral nodes.")

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from NodeListRegistrar import NodeListRegistrar
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

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def generateVoidNodes(generators = "const std::vector<typename %(Dimension)s::Vector>&",
                      Hs = "const std::vector<typename %(Dimension)s::SymTensor>&",
                      mesh = "const Mesh<%(Dimension)s>&",
                      xmin = "const typename %(Dimension)s::Vector&",
                      xmax = "const typename %(Dimension)s::Vector&",
                      numInternal = "const unsigned",
                      nPerh = "const double",
                      threshold = "const double",
                      voidNodes = "NodeList<%(Dimension)s>&"):
    """This algorithm tries to analyze how continuous a node distribution is, and 
if it determines there is an edge to the distribution creates new void nodes
outside that surface.
We assume here that the caller has already created all the boundary ghost 
nodes."""
    return "void"

@PYB11template("Dimension", "moment")
@PYB11implementation("""[](const std::vector<NodeList<%(Dimension)s>*>& nodeLists,
                           const TableKernel<%(Dimension)s>& W,
                           const bool renormalize) -> FieldList<%(Dimension)s, typename MomentTraits<%(Dimension)s, %(moment)s>::Moment> {
                               return nthNodalMoment<%(Dimension)s, typename std::vector<NodeList<%(Dimension)s>*>::const_iterator, %(moment)s>
                               (nodeLists.begin(), nodeLists.end(), W, renormalize);
                           }""")
def nthNodalMoment(nodeLists = "const std::vector<NodeList<%(Dimension)s>*>&",
                   W = "const TableKernel<%(Dimension)s>&",
                   renormalize = "const bool"):
    """ Compute the nth (with n=%(moment)s moment of the local nodal distribution in \\\eta space:

    \\\sum_j  (\\\eta_i)^n W_ij
    -----------------------
         \\\sum_j W_ij
"""
    return "FieldList<%(Dimension)s, typename MomentTraits<%(Dimension)s, %(moment)s>::Moment>"

@PYB11template("Dimension")
@PYB11implementation("""[](const std::vector<NodeList<%(Dimension)s>*>& nodeLists,
                           const TableKernel<%(Dimension)s>& W,
                           const bool useGradientAsKernel,
                           FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>& zerothMoment,
                           FieldList<%(Dimension)s, typename %(Dimension)s::Vector>& firstMoment) {
                                 zerothAndFirstNodalMoments<%(Dimension)s, typename std::vector<NodeList<%(Dimension)s>*>::const_iterator>
                                 (nodeLists.begin(), nodeLists.end(), W, useGradientAsKernel, zerothMoment, firstMoment);
                           }""")
def zerothAndFirstNodalMoments(nodeLists = "const std::vector<NodeList<%(Dimension)s>*>&",
                               W = "const TableKernel<%(Dimension)s>&",
                               useKernelAsGradient = "const bool",
                               zerothMoment = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                               firstMoment = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&"):
    "Compute the non-normalized zeroth and normalized first moment in eta space -- calls nthNodalMoment."
    return "void"

for ndim in dims:
    exec('''
generateVoidNodes%(ndim)id = PYB11TemplateFunction(generateVoidNodes, template_parameters="Dim<%(ndim)i>", pyname="generateVoidNodes")

zerothNodalMoment%(ndim)id = PYB11TemplateFunction(nthNodalMoment, template_parameters=("Dim<%(ndim)i>", "0"), pyname="zerothNodalMoment")
firstNodalMoment%(ndim)id = PYB11TemplateFunction(nthNodalMoment, template_parameters=("Dim<%(ndim)i>", "1"), pyname="firstNodalMoment")
secondNodalMoment%(ndim)id = PYB11TemplateFunction(nthNodalMoment, template_parameters=("Dim<%(ndim)i>", "2"), pyname="secondNodalMoment")

zerothAndFirstNodalMoments%(ndim)id = PYB11TemplateFunction(zerothAndFirstNodalMoments, template_parameters="Dim<%(ndim)i>", pyname="zerothAndFirstNodalMoments")
''' % {"ndim" : ndim})
