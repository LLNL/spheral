"""
Spheral Integrator module.

Provides the base classes and implementations for time integrators in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"DataBase/DataBase.hh"',
            '"DataBase/State.hh"',
            '"DataBase/StateDerivatives.hh"',
            '"Physics/Physics.hh"',
            '"Boundary/Boundary.hh"',
            '"FileIO/FileIO.hh"',
            '"Integrator/Integrator.hh"',
            '"Integrator/PredictorCorrector.hh"',
            '"Integrator/SynchronousRK1.hh"',
            '"Integrator/SynchronousRK2.hh"',
            '"Integrator/SynchronousRK4.hh"',
            '"Integrator/CheapSynchronousRK2.hh"',
            '"Integrator/Verlet.hh"',
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
from Integrator import *
from PredictorCorrectorIntegrator import *

for ndim in dims:
    exec('''
Integrator%(ndim)id = PYB11TemplateClass(Integrator, template_parameters="%(Dimension)s")
PredictorCorrectorIntegrator%(ndim)id = PYB11TemplateClass(PredictorCorrectorIntegrator, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
