"""
Spheral SPH module.

Provides implementations of SPH, PSPH, and ASPH 
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

from NodeCoupling import *
from SPHHydroBase import *
from PSPHHydroBase import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"SPH/SPHHydroBase.hh"',
            '"SPH/PSPHHydroBase.hh"',
            '"SPH/computeSPHSumMassDensity.hh"',
            '"SPH/computeSPHOmegaGradhCorrection.hh"',
            '"SPH/SPHHydroBaseRZ.hh"',
            '"SPH/SPHHydroBaseGSRZ.hh"',
            '"SPH/SolidSPHHydroBase.hh"',
            '"SPH/SolidSPHHydroBaseRZ.hh"',
            '"SPH/NodeCoupling.hh"',
            '"SPH/DamagedNodeCoupling.hh"',
            '"SPH/DamagedNodeCouplingWithFrags.hh"',
            '"FileIO/FileIO.hh"',
            '"ArtificialViscosity/ArtificialViscosity.hh"',
            '<vector>',
            '<string>',
            '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
DamagedNodeCoupling%(ndim)id = PYB11TemplateClass(DamagedNodeCoupling, template_parameters="%(Dimension)s")
DamagedNodeCouplingWithFrags%(ndim)id = PYB11TemplateClass(DamagedNodeCouplingWithFrags, template_parameters="%(Dimension)s")

SPHHydroBase%(ndim)id = PYB11TemplateClass(SPHHydroBase, template_parameters="%(Dimension)s")
PSPHHydroBase%(ndim)id = PYB11TemplateClass(PSPHHydroBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
