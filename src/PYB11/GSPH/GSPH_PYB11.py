"""
Spheral GSPH module.

Provides implementations of Riemann SPH
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from GenericRiemannHydro import *
from GSPHHydroBase import *
from MFMHydroBase import *
from WaveSpeeds import *
from Limiters import *
from RiemannSolvers import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"GSPH/GenericRiemannHydro.hh"',
                  '"GSPH/GSPHHydroBase.hh"',
                  '"GSPH/MFMHydroBase.hh"',
                  '"GSPH/WaveSpeeds/WaveSpeedBase.hh"',
                  '"GSPH/WaveSpeeds/AcousticWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/DavisWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/EinfeldtWaveSpeed.hh"',
                  '"GSPH/Limiters/LimiterBase.hh"',
                  '"GSPH/Limiters/MinModLimiter.hh"',
                  '"GSPH/Limiters/VanLeerLimiter.hh"',
                  '"GSPH/Limiters/VanAlbaLimiter.hh"',
                  '"GSPH/Limiters/SuperbeeLimiter.hh"',
                  '"GSPH/Limiters/OspreLimiter.hh"',
                  '"GSPH/RiemannSolvers/RiemannSolverBase.hh"',
                  '"GSPH/RiemannSolvers/HLLC.hh"',
                  '"GSPH/RiemannSolvers/SecondOrderArtificialViscosity.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
GradientType = PYB11enum(("RiemannGradient",
                          "HydroAccelerationGradient",
                          "SPHGradient",
                          "MixedMethodGradient",
                          "SPHSameTimeGradient",
                          "SPHUncorrectedGradient",
                          "NoGradient"), export_values = True)

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
GenericRiemannHydro%(ndim)id = PYB11TemplateClass(GenericRiemannHydro, template_parameters="%(Dimension)s")
GSPHHydroBase%(ndim)id = PYB11TemplateClass(GSPHHydroBase, template_parameters="%(Dimension)s")
MFMHydroBase%(ndim)id = PYB11TemplateClass(MFMHydroBase, template_parameters="%(Dimension)s")
WaveSpeedBase%(ndim)id = PYB11TemplateClass(WaveSpeedBase, template_parameters="%(Dimension)s")
AcousticWaveSpeed%(ndim)id = PYB11TemplateClass(AcousticWaveSpeed, template_parameters="%(Dimension)s")
DavisWaveSpeed%(ndim)id = PYB11TemplateClass(DavisWaveSpeed, template_parameters="%(Dimension)s")
EinfeldtWaveSpeed%(ndim)id = PYB11TemplateClass(EinfeldtWaveSpeed, template_parameters="%(Dimension)s")
LimiterBase%(ndim)id = PYB11TemplateClass(LimiterBase, template_parameters="%(Dimension)s")
MinModLimiter%(ndim)id = PYB11TemplateClass(MinModLimiter, template_parameters="%(Dimension)s")
VanLeerLimiter%(ndim)id = PYB11TemplateClass(VanLeerLimiter, template_parameters="%(Dimension)s")
VanAlbaLimiter%(ndim)id = PYB11TemplateClass(VanAlbaLimiter, template_parameters="%(Dimension)s")
SuperbeeLimiter%(ndim)id = PYB11TemplateClass(SuperbeeLimiter, template_parameters="%(Dimension)s")
OspreLimiter%(ndim)id = PYB11TemplateClass(OspreLimiter, template_parameters="%(Dimension)s")
RiemannSolverBase%(ndim)id = PYB11TemplateClass(RiemannSolverBase, template_parameters="%(Dimension)s")
HLLC%(ndim)id = PYB11TemplateClass(HLLC, template_parameters="%(Dimension)s")
SecondOrderArtificialViscosity%(ndim)id = PYB11TemplateClass(SecondOrderArtificialViscosity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

