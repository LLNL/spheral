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
PYB11includes += ['"ArtificialViscosity/ArtificialViscosityHandle.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
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
from ArtificialViscosityHandle import *
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
    Dimension = f"Dim<{ndim}>"
    exec(f'''
ArtificialViscosityHandle{ndim}d = PYB11TemplateClass(ArtificialViscosityHandle, template_parameters="{Dimension}")
ScalarArtificialViscosity{ndim}d = PYB11TemplateClass(ArtificialViscosity, template_parameters=("{Dimension}", "{Dimension}::Scalar"))
TensorArtificialViscosity{ndim}d = PYB11TemplateClass(ArtificialViscosity, template_parameters=("{Dimension}", "{Dimension}::Tensor"))
MonaghanGingoldViscosity{ndim}d = PYB11TemplateClass(MonaghanGingoldViscosity, template_parameters="{Dimension}")
TensorMonaghanGingoldViscosity{ndim}d = PYB11TemplateClass(TensorMonaghanGingoldViscosity, template_parameters="{Dimension}")
LimitedMonaghanGingoldViscosity{ndim}d = PYB11TemplateClass(LimitedMonaghanGingoldViscosity, template_parameters="{Dimension}")
MorrisMonaghanReducingViscosity{ndim}d = PYB11TemplateClass(MorrisMonaghanReducingViscosity, template_parameters="{Dimension}")
CullenDehnenViscosity{ndim}d = PYB11TemplateClass(CullenDehnenViscosity, template_parameters="{Dimension}")
FiniteVolumeViscosity{ndim}d = PYB11TemplateClass(FiniteVolumeViscosity, template_parameters="{Dimension}")
TensorSVPHViscosity{ndim}d = PYB11TemplateClass(TensorSVPHViscosity, template_parameters="{Dimension}")
TensorCRKSPHViscosity{ndim}d = PYB11TemplateClass(TensorCRKSPHViscosity, template_parameters="{Dimension}")
''')
