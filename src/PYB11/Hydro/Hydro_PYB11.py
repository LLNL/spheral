"""
Spheral Hydro module.

Provides the support classes for hydro algorithms.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"DataBase/DataBase.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Hydro/HydroFieldNames.hh"',
                  '"Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"',
                  '"Hydro/EntropyPolicy.hh"',
                  '"Hydro/GammaPolicy.hh"',
                  '"Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"',
                  '"Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"',
                  '"Hydro/SphericalPositionPolicy.hh"',
                  '"Hydro/PressurePolicy.hh"',
                  '"Hydro/SoundSpeedPolicy.hh"',
                  '"Hydro/SpecificFromTotalThermalEnergyPolicy.hh"',
                  '"Hydro/SpecificThermalEnergyPolicy.hh"',
                  '"Hydro/SumVoronoiMassDensityPolicy.hh"',
                  '"Hydro/VolumePolicy.hh"',
                  '"Hydro/VoronoiMassDensityPolicy.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from HydroFieldNames import *
from CompatibleDifferenceSpecificThermalEnergyPolicy import *
from NonSymmetricSpecificThermalEnergyPolicy import *
from SpecificThermalEnergyPolicy import *
from SpecificFromTotalThermalEnergyPolicy import *
from EntropyPolicy import *
from GammaPolicy import *
from PressurePolicy import *
from SoundSpeedPolicy import *
from SumVoronoiMassDensityPolicy import *
from VoronoiMassDensityPolicy import *
from VolumePolicy import *

for ndim in dims:

    exec(f'''
CompatibleDifferenceSpecificThermalEnergyPolicy{ndim}d = PYB11TemplateClass(CompatibleDifferenceSpecificThermalEnergyPolicy, template_parameters="Dim<{ndim}>")
NonSymmetricSpecificThermalEnergyPolicy{ndim}d = PYB11TemplateClass(NonSymmetricSpecificThermalEnergyPolicy, template_parameters="Dim<{ndim}>")
SpecificThermalEnergyPolicy{ndim}d = PYB11TemplateClass(SpecificThermalEnergyPolicy, template_parameters="Dim<{ndim}>")
SpecificFromTotalThermalEnergyPolicy{ndim}d = PYB11TemplateClass(SpecificFromTotalThermalEnergyPolicy, template_parameters="Dim<{ndim}>")
EntropyPolicy{ndim}d = PYB11TemplateClass(EntropyPolicy, template_parameters="Dim<{ndim}>")
GammaPolicy{ndim}d = PYB11TemplateClass(GammaPolicy, template_parameters="Dim<{ndim}>")
PressurePolicy{ndim}d = PYB11TemplateClass(PressurePolicy, template_parameters="Dim<{ndim}>")
SoundSpeedPolicy{ndim}d = PYB11TemplateClass(SoundSpeedPolicy, template_parameters="Dim<{ndim}>")
SumVoronoiMassDensityPolicy{ndim}d = PYB11TemplateClass(SumVoronoiMassDensityPolicy, template_parameters="Dim<{ndim}>")
VoronoiMassDensityPolicy{ndim}d = PYB11TemplateClass(VoronoiMassDensityPolicy, template_parameters="Dim<{ndim}>")
VolumePolicy{ndim}d = PYB11TemplateClass(VolumePolicy, template_parameters="Dim<{ndim}>")
''')

#-------------------------------------------------------------------------------
# Specialized types
#-------------------------------------------------------------------------------
if 1 in dims:
    from SphericalPositionPolicy import *

if 2 in dims:
    from RZNonSymmetricSpecificThermalEnergyPolicy import *
