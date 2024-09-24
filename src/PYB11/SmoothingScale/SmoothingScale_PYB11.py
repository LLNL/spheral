"""
Spheral SmoothignScale module.

Provides the physics packages for updating the smoothing scale
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SmoothingScale/SmoothingScaleBase.hh"',
                  '"SmoothingScale/FixedSmoothingScale.hh"',
                  '"SmoothingScale/SPHSmoothingScale.hh"',
                  '"SmoothingScale/ASPHSmoothingScale.hh"',
                  '"SmoothingScale/ASPHSmoothingScaleUserFilter.hh"',
                  '"Kernel/TableKernel.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
HEvolutionType = PYB11enum(("IdealH", "IntegrateH", "FixedH"), export_values=True,
                     doc="The choices for updating the smoothing scale.")

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from SmoothingScaleBase import SmoothingScaleBase
from FixedSmoothingScale import FixedSmoothingScale
from SPHSmoothingScale import SPHSmoothingScale
from ASPHSmoothingScale import ASPHSmoothingScale
from ASPHSmoothingScaleUserFilter import ASPHSmoothingScaleUserFilter

for ndim in dims:
    exec(f'''
SmoothingScaleBase{ndim}d = PYB11TemplateClass(SmoothingScaleBase, template_parameters="Dim<{ndim}>")
FixedSmoothingScale{ndim}d = PYB11TemplateClass(FixedSmoothingScale, template_parameters="Dim<{ndim}>")
SPHSmoothingScale{ndim}d = PYB11TemplateClass(SPHSmoothingScale, template_parameters="Dim<{ndim}>")
ASPHSmoothingScale{ndim}d = PYB11TemplateClass(ASPHSmoothingScale, template_parameters="Dim<{ndim}>")
ASPHSmoothingScaleUserFilter{ndim}d = PYB11TemplateClass(ASPHSmoothingScaleUserFilter, template_parameters="Dim<{ndim}>")
''')
