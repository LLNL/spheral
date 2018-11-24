"""
Spheral Integrator module.

Provides the base classes and implementations for time integrators in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DataBase/DataBase.hh"',
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
                  '"Integrator/Verlet.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from Integrator import *
from PredictorCorrectorIntegrator import *
from SynchronousRK1Integrator import *
from SynchronousRK2Integrator import *
from SynchronousRK4Integrator import *
from CheapSynchronousRK2Integrator import *
from VerletIntegrator import *

for ndim in dims:
    exec('''
Integrator%(ndim)id = PYB11TemplateClass(Integrator, template_parameters="%(Dimension)s")
PredictorCorrectorIntegrator%(ndim)id = PYB11TemplateClass(PredictorCorrectorIntegrator, template_parameters="%(Dimension)s")
SynchronousRK1Integrator%(ndim)id = PYB11TemplateClass(SynchronousRK1Integrator, template_parameters="%(Dimension)s")
SynchronousRK2Integrator%(ndim)id = PYB11TemplateClass(SynchronousRK2Integrator, template_parameters="%(Dimension)s")
SynchronousRK4Integrator%(ndim)id = PYB11TemplateClass(SynchronousRK4Integrator, template_parameters="%(Dimension)s")
CheapSynchronousRK2Integrator%(ndim)id = PYB11TemplateClass(CheapSynchronousRK2Integrator, template_parameters="%(Dimension)s")
VerletIntegrator%(ndim)id = PYB11TemplateClass(VerletIntegrator, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
