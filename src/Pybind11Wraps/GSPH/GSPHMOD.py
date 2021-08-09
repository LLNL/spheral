"""
Spheral GSPH module.

Provides implementations of Riemann SPH
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from GSPHHydroBase import *
from WaveSpeeds import *
from Limiters import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"GSPH/GSPHHydroBase.hh"',
                  '"GSPH/WaveSpeeds/WaveSpeedBase.hh"',
                  '"GSPH/WaveSpeeds/AcousticWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/DavisWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/EinfeldtWaveSpeed.hh"',
                  '"GSPH/Limiters/SlopeLimiterBase.hh"',
                  '"GSPH/Limiters/MinModLimiter.hh"',
                  '"GSPH/Limiters/VanLeerLimiter.hh"',
                  '"GSPH/Limiters/VanAlbaLimiter.hh"',
                  '"GSPH/Limiters/SuperbeeLimiter.hh"',
                  '"GSPH/Limiters/OspreLimiter.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
GSPHHydroBase%(ndim)id = PYB11TemplateClass(GSPHHydroBase, template_parameters="%(Dimension)s")
WaveSpeedBase%(ndim)id = PYB11TemplateClass(WaveSpeedBase, template_parameters="%(Dimension)s")
AcousticWaveSpeed%(ndim)id = PYB11TemplateClass(AcousticWaveSpeed, template_parameters="%(Dimension)s")
DavisWaveSpeed%(ndim)id = PYB11TemplateClass(DavisWaveSpeed, template_parameters="%(Dimension)s")
EinfeldtWaveSpeed%(ndim)id = PYB11TemplateClass(EinfeldtWaveSpeed, template_parameters="%(Dimension)s")
SlopeLimiterBase%(ndim)id = PYB11TemplateClass(SlopeLimiterBase, template_parameters="%(Dimension)s")
MinModLimiter%(ndim)id = PYB11TemplateClass(MinModLimiter, template_parameters="%(Dimension)s")
VanLeerLimiter%(ndim)id = PYB11TemplateClass(VanLeerLimiter, template_parameters="%(Dimension)s")
VanAlbaLimiter%(ndim)id = PYB11TemplateClass(VanAlbaLimiter, template_parameters="%(Dimension)s")
SuperbeeLimiter%(ndim)id = PYB11TemplateClass(SuperbeeLimiter, template_parameters="%(Dimension)s")
OspreLimiter%(ndim)id = PYB11TemplateClass(OspreLimiter, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

