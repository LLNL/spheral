"""
Spheral ArtificialViscosity module.

Provides the artificial viscosity algorithms for use with the hydrodynamics methods.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"ArtificialViscosity/ArtificialViscosity.hh"',
                  '"ArtificialViscosity/MonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/MorrisMonaghanReducingViscosity.hh"',
                  '"ArtificialViscosity/CullenDehnenViscosity.hh"',
                  '"ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/FiniteVolumeViscosity.hh"',
                  '"ArtificialViscosity/TensorSVPHViscosity.hh"',
                  '"ArtificialViscosity/TensorCRKSPHViscosity.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from ArtificialViscosity import *
from MonaghanGingoldViscosity import *
from TensorMonaghanGingoldViscosity import *
from LimitedMonaghanGingoldViscosity import *
from MorrisMonaghanReducingViscosity import *
from CullenDehnenViscosity import *
from FiniteVolumeViscosity import *
from TensorSVPHViscosity import *
from TensorCRKSPHViscosity import *

for ndim in dims:
    exec('''
ArtificialViscosity%(ndim)id = PYB11TemplateClass(ArtificialViscosity, template_parameters="%(Dimension)s")
MonaghanGingoldViscosity%(ndim)id = PYB11TemplateClass(MonaghanGingoldViscosity, template_parameters="%(Dimension)s")
TensorMonaghanGingoldViscosity%(ndim)id = PYB11TemplateClass(TensorMonaghanGingoldViscosity, template_parameters="%(Dimension)s")
LimitedMonaghanGingoldViscosity%(ndim)id = PYB11TemplateClass(LimitedMonaghanGingoldViscosity, template_parameters="%(Dimension)s")
MorrisMonaghanReducingViscosity%(ndim)id = PYB11TemplateClass(MorrisMonaghanReducingViscosity, template_parameters="%(Dimension)s")
CullenDehnenViscosity%(ndim)id = PYB11TemplateClass(CullenDehnenViscosity, template_parameters="%(Dimension)s")
FiniteVolumeViscosity%(ndim)id = PYB11TemplateClass(FiniteVolumeViscosity, template_parameters="%(Dimension)s")
TensorSVPHViscosity%(ndim)id = PYB11TemplateClass(TensorSVPHViscosity, template_parameters="%(Dimension)s")
TensorCRKSPHViscosity%(ndim)id = PYB11TemplateClass(TensorCRKSPHViscosity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})
