"""
Spheral Physics module.

Provides the Physics, GenericHydro, and other physics base classes
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Physics/Physics.hh"',
            '"Physics/GenericHydro.hh"',
            '"Physics/GenericBodyForce.hh"',
            '"Boundary/Boundary.hh"',
            '"ArtificialViscosity/ArtificialViscosity.hh"',
            '"Kernel/TableKernel.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
MassDensityType = PYB11enum(("SumDensity", 
                             "RigorousSumDensity",
                             "HybridSumDensity",
                             "IntegrateDensity",
                             "VoronoiCellDensity",
                             "SumVoronoiCellDensity",
                             "CorrectedSumDensity"), export_values=True)
HEvolutionType = PYB11enum(("IdealH", 
                            "IntegrateH"), export_values = True)

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from Physics import *
from GenericHydro import *

for ndim in dims:
    exec('''
Physics%(ndim)id = PYB11TemplateClass(Physics, template_parameters="%(Dimension)s")
GenericHydro%(ndim)id = PYB11TemplateClass(GenericHydro, template_parameters="%(Dimension)s")

vector_of_Physics%(ndim)id = PYB11_bind_vector("Physics<%(Dimension)s>*", opaque=True)
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})
