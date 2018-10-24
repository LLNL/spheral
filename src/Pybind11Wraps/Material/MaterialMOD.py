"""
Spheral Material module.

Provides the interfaces and basic fluid material models.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Material/PhysicalConstants.hh"',
            '"Material/EquationOfState.hh"',
            '"Material/GammaLawGas.hh"',
            '"Material/PolytropicEquationOfState.hh"',
            '"Material/IsothermalEquationOfState.hh"',
            '"Field/Field.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
MaterialPressureMinType = PYB11enum(("PressureFloor", "ZeroPressure"), export_values=True,
                                    doc="Choices for how to apply minimum pressure threshold")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from PhysicalConstants import *
from EquationOfState import *
from GammaLawGas import *
from PolytropicEquationOfState import *
from IsothermalEquationOfState import *

for ndim in dims:
    exec('''
EquationOfState%(ndim)id = PYB11TemplateClass(EquationOfState, template_parameters="Dim<%(ndim)i>")
GammaLawGas%(ndim)id = PYB11TemplateClass(GammaLawGas, template_parameters="Dim<%(ndim)i>")
PolytropicEquationOfState%(ndim)id = PYB11TemplateClass(PolytropicEquationOfState, template_parameters="Dim<%(ndim)i>")
IsothermalEquationOfState%(ndim)id = PYB11TemplateClass(IsothermalEquationOfState, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})
