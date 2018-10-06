"""
Spheral FieldOperations module.

Functions for assorted operation over Fields/FieldLists (interpolation, 
resampling, gradient, etc.)
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Field/Field.hh"',
            '"Field/FieldList.hh"',
            '"Kernel/TableKernel.hh"',
            '"FieldOperations/FieldListFunctions.hh"',
            '"FieldOperations/FieldListFunctionsMash.hh"',
            '"FieldOperations/FieldListSecondDerivatives.hh"',
            '"FieldOperations/PairWiseFieldListFunctions.hh"',
            '"FieldOperations/sampleMultipleFields2Lattice.hh"',
            '"FieldOperations/binFieldList2Lattice.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Define the templated methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
def smoothFields(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                 position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                 weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                 mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                 rho = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                 Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                 kernel = "const TableKernel<%(Dimension)s>&"):
    "Calculate a smoothed estimate of the given FieldList."
    return "FieldList<%(Dimension)s, %(DataType)s>"

@PYB11template("Dimension", "DataType")
def gradient(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
             position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
             weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
             mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
             rho = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
             Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
             kernel = "const TableKernel<%(Dimension)s>&"):
    "Calculate the gradient of the given FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>"

@PYB11template("Dimension", "DataType")
def divergence(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
               position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
               weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
               mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
               rho = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
               Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
               kernel = "const TableKernel<%(Dimension)s>&"):
    "Calculate the divergence of the given FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::DivergenceType>"

@PYB11template("Dimension", "DataType")
def limiter(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
            gradient = "const FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>&",
            position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
            Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
            kernel = "const TableKernel<%(Dimension)s>&"):
    "Calculate a tensor limiter appropriate for a monotonically limited gradient."
    return "FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>"

#-------------------------------------------------------------------------------
# Instantiate the per dimension bindings
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
smoothScalarFields%(ndim)id = PYB11TemplateFunction(smoothFields, template_parameters=("%(Dimension)s", "%(Scalar)s"))
smoothVectorFields%(ndim)id = PYB11TemplateFunction(smoothFields, template_parameters=("%(Dimension)s", "%(Vector)s"))
smoothTensorFields%(ndim)id = PYB11TemplateFunction(smoothFields, template_parameters=("%(Dimension)s", "%(Tensor)s"))
smoothSymTensorFields%(ndim)id = PYB11TemplateFunction(smoothFields, template_parameters=("%(Dimension)s", "%(SymTensor)s"))

gradientScalar%(ndim)id = PYB11TemplateFunction(gradient, template_parameters=("%(Dimension)s", "%(Scalar)s"), pyname="gradient")
gradientVector%(ndim)id = PYB11TemplateFunction(gradient, template_parameters=("%(Dimension)s", "%(Vector)s"), pyname="gradient")

divergenceVector%(ndim)id = PYB11TemplateFunction(divergence, template_parameters=("%(Dimension)s", "%(Vector)s"), pyname="divergence")
divergenceTensor%(ndim)id = PYB11TemplateFunction(divergence, template_parameters=("%(Dimension)s", "%(Tensor)s"), pyname="divergence")
divergenceSymTensor%(ndim)id = PYB11TemplateFunction(divergence, template_parameters=("%(Dimension)s", "%(SymTensor)s"), pyname="divergence")

scalarLimiter%(ndim)id = PYB11TemplateFunction(limiter, template_parameters=("%(Dimension)s", "%(Scalar)s"), pyname="scalarLimiter%(ndim)id")
vectorLimiter%(ndim)id = PYB11TemplateFunction(limiter, template_parameters=("%(Dimension)s", "%(Vector)s"), pyname="vectorLimiter%(ndim)id")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "Scalar"    : "Dim<" + str(ndim) + ">::Scalar",
       "Vector"    : "Dim<" + str(ndim) + ">::Vector",
       "Tensor"    : "Dim<" + str(ndim) + ">::Tensor",
       "SymTensor" : "Dim<" + str(ndim) + ">::SymTensor"})
