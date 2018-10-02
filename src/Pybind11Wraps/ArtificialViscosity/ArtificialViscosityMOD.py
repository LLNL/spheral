"""
Spheral ArtificialViscosity module.

Provides the artificial viscosity algorithms for use with the hydrodynamics methods.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"ArtificialViscosity/ArtificialViscosity.hh"',
            '"ArtificialViscosity/MonaghanGingoldViscosity.hh"',
            '"ArtificialViscosity/CRKSPHMonaghanGingoldViscosity.hh"',
            '"ArtificialViscosity/MorrisMonaghanReducingViscosity.hh"',
            '"ArtificialViscosity/CullenDehnenViscosity.hh"',
            '"ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"',
            '"ArtificialViscosity/FiniteVolumeViscosity.hh"',
            '"ArtificialViscosity/TensorSVPHViscosity.hh"',
            '"ArtificialViscosity/TensorCRKSPHViscosity.hh"',
            '"ArtificialViscosity/VonNeumanViscosity.hh"',
            '"ArtificialViscosity/MonaghanGingoldViscosityGSRZ.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from ArtificialViscosity import *

for ndim in dims:
    exec('''
ArtificialViscosity%(ndim)id = PYB11TemplateClass(ArtificialViscosity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})
