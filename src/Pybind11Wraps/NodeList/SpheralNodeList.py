"""
Spheral NodeList module.

Provides the fundamental NodeList classes.
"""

from PYB11Generator import *
from spheralDimensions import spheralDimensions as PYB11dimensions
dims = PYB11dimensions()

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
            '"Mesh/Mesh.hh"']

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
# NodeLists
#-------------------------------------------------------------------------------
from NodeList import NodeList

NodeList1d = PYB11TemplateClass(NodeList,
                                template_parameters = {"Dimension" : "Dim<1>",
                                                       "Scalar"    : "Dim<1>::Scalar",
                                                       "Vector"    : "Dim<1>::Vector",
                                                       "Tensor"    : "Dim<1>::Tensor",
                                                       "SymTensor" : "Dim<1>::SymTensor"})

# for ndim in dims:
#     exec('''
# NodeList%(ndim)id =  PYB11TemplateClass(NodeList, pyname="NodeList%(ndim)id",
#                                         template_parameters=("Dim<%(ndim)i>"))
# ''' % {"ndim" : ndim,
#        "Scalar" : "Dim<%i>::Scalar" % ndim,
#        "Vector" : "Dim<%i>::Vector" % ndim,
#        "Tensor" : "Dim<%i>::Tensor" % ndim,
#        "SymTensor" : "Dim<%i>::Symensor" % ndim})
