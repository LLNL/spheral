"""
Spheral KernelIntegrator module
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from BilinearIndex import *
from FlatConnectivity import *
from KernelIntegrationData import *
from IntegrationCoefficient import *
from IntegrationKernel import *
from KernelIntegral import *
from KernelIntegrator import *
from ManufacturedSolution import *
from RKIntegrationKernel import *
from SPHIntegrationKernel import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['<memory>',
                  '<string>',
                  '<tuple>',
                  '<vector>',
                  '<unordered_map>',
                  '"Boundary/Boundary.hh"',
                  '"DataBase/DataBase.hh"',
                  '"DataBase/State.hh"',
                  '"Field/FieldList.hh"',
                  '"Geometry/GeomPlane.hh"',
                  '"Hydro/HydroFieldNames.hh"',
                  '"Registrar/BoundaryDependent.hh"',
                  '"Registrar/BoundaryRegistrar.hh"',
                  '"Registrar/CombinedRegistrar.hh"',
                  '"Registrar/ConnectivityDependent.hh"',
                  '"Registrar/ConnectivityRegistrar.hh"',
                  '"Registrar/StateDependent.hh"',
                  '"Registrar/StateRegistrar.hh"',
                  '"Registrar/UpdateDependent.hh"',
                  '"Registrar/UpdateRegistrar.hh"',
                  '"KernelIntegrator/BilinearIndex.hh"',
                  '"KernelIntegrator/FlatConnectivity.hh"',
                  '"KernelIntegrator/IntegrationCoefficient.hh"',
                  '"KernelIntegrator/IntegrationKernel.hh"',
                  '"KernelIntegrator/KernelIntegral.hh"',
                  '"KernelIntegrator/KernelIntegrationData.hh"',
                  '"KernelIntegrator/KernelIntegrator.hh"',
                  '"KernelIntegrator/ManufacturedSolution.hh"',
                  '"KernelIntegrator/RKIntegrationKernel.hh"',
                  '"KernelIntegrator/SPHIntegrationKernel.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiations
#-------------------------------------------------------------------------------
for ndim in dims:
    # Dimension-dependent
    exec('''
BilinearIndex%(ndim)id = PYB11TemplateClass(BilinearIndex, template_parameters="Dim<%(ndim)i>")
FlatConnectivity%(ndim)id = PYB11TemplateClass(FlatConnectivity, template_parameters="Dim<%(ndim)i>")
IntegrationKernel%(ndim)id = PYB11TemplateClass(IntegrationKernel, template_parameters="Dim<%(ndim)i>")
SPHIntegrationKernel%(ndim)id = PYB11TemplateClass(SPHIntegrationKernel, template_parameters="Dim<%(ndim)i>")
KernelIntegrator%(ndim)id = PYB11TemplateClass(KernelIntegrator, template_parameters="Dim<%(ndim)i>")
KernelIntegrationData%(ndim)id = PYB11TemplateClass(KernelIntegrationData, template_parameters="Dim<%(ndim)i>")
KernelIntegralBase%(ndim)id = PYB11TemplateClass(KernelIntegralBase, template_parameters="Dim<%(ndim)i>")
LinearKernel%(ndim)id = PYB11TemplateClass(LinearKernel, template_parameters="Dim<%(ndim)i>")
LinearGrad%(ndim)id = PYB11TemplateClass(LinearGrad, template_parameters="Dim<%(ndim)i>")
LinearKernelVector%(ndim)id = PYB11TemplateClass(LinearKernelVector, template_parameters="Dim<%(ndim)i>")
LinearKernelStdVector%(ndim)id = PYB11TemplateClass(LinearKernelStdVector, template_parameters="Dim<%(ndim)i>")
LinearGradStdVector%(ndim)id = PYB11TemplateClass(LinearGradStdVector, template_parameters="Dim<%(ndim)i>")
BilinearKernelKernel%(ndim)id = PYB11TemplateClass(BilinearKernelKernel, template_parameters="Dim<%(ndim)i>")
BilinearGradKernel%(ndim)id = PYB11TemplateClass(BilinearGradKernel, template_parameters="Dim<%(ndim)i>")
BilinearKernelGrad%(ndim)id = PYB11TemplateClass(BilinearKernelGrad, template_parameters="Dim<%(ndim)i>")
BilinearGradDotGrad%(ndim)id = PYB11TemplateClass(BilinearGradDotGrad, template_parameters="Dim<%(ndim)i>")
BilinearGradProdGrad%(ndim)id = PYB11TemplateClass(BilinearGradProdGrad, template_parameters="Dim<%(ndim)i>")
BilinearSurfaceNormalKernelKernelFromGrad%(ndim)id = PYB11TemplateClass(BilinearSurfaceNormalKernelKernelFromGrad, template_parameters="Dim<%(ndim)i>")
LinearSurfaceNormalKernel%(ndim)id = PYB11TemplateClass(LinearSurfaceNormalKernel, template_parameters="Dim<%(ndim)i>") 
LinearSurfaceKernel%(ndim)id = PYB11TemplateClass(LinearSurfaceKernel, template_parameters="Dim<%(ndim)i>")
LinearSurfaceNormalKernelStdVector%(ndim)id = PYB11TemplateClass(LinearSurfaceNormalKernelStdVector, template_parameters="Dim<%(ndim)i>") 
BilinearSurfaceKernelKernel%(ndim)id = PYB11TemplateClass(BilinearSurfaceKernelKernel, template_parameters="Dim<%(ndim)i>") 
BilinearSurfaceNormalKernelKernel%(ndim)id = PYB11TemplateClass(BilinearSurfaceNormalKernelKernel, template_parameters="Dim<%(ndim)i>") 
BilinearSurfaceNormalKernelDotGrad%(ndim)id = PYB11TemplateClass(BilinearSurfaceNormalKernelDotGrad, template_parameters="Dim<%(ndim)i>") 
CellCoefficient%(ndim)id = PYB11TemplateClass(CellCoefficient, template_parameters="Dim<%(ndim)i>")
SurfaceNormalCoefficient%(ndim)id = PYB11TemplateClass(SurfaceNormalCoefficient, template_parameters="Dim<%(ndim)i>") 
ManufacturedFunction%(ndim)id = PYB11TemplateClass(ManufacturedFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedSteadyStateFunction%(ndim)id = PYB11TemplateClass(ManufacturedSteadyStateFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedConstantFunction%(ndim)id = PYB11TemplateClass(ManufacturedConstantFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedSinusoidalFunction%(ndim)id = PYB11TemplateClass(ManufacturedSinusoidalFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedWaveFunction%(ndim)id = PYB11TemplateClass(ManufacturedWaveFunction, template_parameters="Dim<%(ndim)i>")
ManufacturedTransportSolution%(ndim)id = PYB11TemplateClass(ManufacturedTransportSolution, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    # Dependent on primitives
    dim_types = (("Dim<%i>::Scalar" % ndim, "Scalar"),
                 ("Dim<%i>::Vector" % ndim, "Vector"),
                 ("Dim<%i>::SymTensor" % ndim, "SymTensor"),
                 ("Dim<%i>::Tensor" % ndim, "Tensor"),
                 ("std::vector<Dim<%i>::Scalar>" % ndim, "StdVectorScalar"),
                 ("std::vector<Dim<%i>::Vector>" % ndim, "StdVectorVector"))
    for (value, label) in dim_types:
        exec('''
%(label)sConstantIntegrationCoefficient%(ndim)id = PYB11TemplateClass(ConstantIntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sIntegrationCoefficient%(ndim)id = PYB11TemplateClass(IntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sFieldListIntegrationCoefficient%(ndim)id = PYB11TemplateClass(FieldListIntegrationCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sIntegralDependsOnCoefficient%(ndim)id = PYB11TemplateClass(IntegralDependsOnCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sIntegralDependsOnFieldListCoefficient%(ndim)id = PYB11TemplateClass(IntegralDependsOnFieldListCoefficient, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sLinearIntegral%(ndim)id = PYB11TemplateClass(LinearIntegral, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sBilinearIntegral%(ndim)id = PYB11TemplateClass(BilinearIntegral, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sLinearSurfaceDependentIntegral%(ndim)id = PYB11TemplateClass(LinearSurfaceDependentIntegral, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
%(label)sBilinearSurfaceDependentIntegral%(ndim)id = PYB11TemplateClass(BilinearSurfaceDependentIntegral, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
        if not "StdVector" in label:
            exec('''
%(label)sBilinearMultiplyByFieldList%(ndim)id = PYB11TemplateClass(BilinearMultiplyByFieldList, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
    
    # Dependent on data types
    data_types = dim_types
    for (value, label) in dim_types:
        data_types += (("std::vector<%(value)s>" % {"value" : value},
                        "StdVector%(label)s" % {"label" : label}),)
        
        
    for (value, label) in data_types:
        exec('''
%(label)sKernelIntegral%(ndim)id = PYB11TemplateClass(KernelIntegral, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    # Dependent on order
    for order in (0, 1, 2, 3, 4, 5, 6, 7):
        exec('''
RKIntegrationKernel%(ndim)id%(order)s = PYB11TemplateClass(RKIntegrationKernel, template_parameters=("Dim<%(ndim)i>", "%(order)s"))
''' % {"ndim" : ndim,
       "order" : order})
