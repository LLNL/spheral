"""
Spheral SVPH module.

Provides implementations of SVPH (Smoothed Voronoi Particle Hydrodynamics)
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from SVPHFieldNames import *
from SVPHFacetedHydroBase import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SVPH/SVPHFieldNames.hh"',
                  '"SVPH/sampleFieldListSVPH.hh"',
                  '"SVPH/gradientFieldListSVPH.hh"',
                  '"SVPH/SVPHHydroBase.hh"',
                  '"SVPH/SVPHFacetedHydroBase.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
def sampleFieldListSVPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                        position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                        Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                        connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                        W = "const TableKernel<%(Dimension)s>&",
                        mesh = "const Mesh<%(Dimension)s>&",
                        firstOrderConsistent = "const bool"):
    "Use SVPH to sample a FieldList."
    return "FieldList<%(Dimension)s, %(DataType)s>"

@PYB11template("Dimension", "DataType")
def gradientFieldListSVPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                          position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                          Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                          connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                          W = "const TableKernel<%(Dimension)s>&",
                          mesh = "const Mesh<%(Dimension)s>&",
                          firstOrderConsistent = "const bool"):
    "Use SVPH to take the gradient of a FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
SVPHFacetedHydroBase%(ndim)id = PYB11TemplateClass(SVPHFacetedHydroBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

    # SVPH interpolation
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim,
                    "Dim<%i>::Tensor" % ndim,
                    "Dim<%i>::SymTensor" % ndim):
        exec('''
sampleFieldListSVPH%(label)s = PYB11TemplateFunction(sampleFieldListSVPH, template_parameters=("%(Dimension)s", "%(element)s"), pyname="sampleFieldListSVPH")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

    # SVPH gradient
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim):
        exec('''
gradientFieldListSVPH%(label)s = PYB11TemplateFunction(gradientFieldListSVPH, template_parameters=("%(Dimension)s", "%(element)s"), pyname="gradientFieldListSVPH")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

