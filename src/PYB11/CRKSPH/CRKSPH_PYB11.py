"""
Spheral CRKSPH module.

Provides implementations of CRKSPH (fluid and solid).
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from CRKSPHBase import *
from CRKSPH import *
from SolidCRKSPH import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"CRKSPH/CRKSPHBase.hh"',
                  '"CRKSPH/CRKSPH.hh"',
                  '"CRKSPH/CRKSPHRZ.hh"',
                  '"CRKSPH/SolidCRKSPH.hh"',
                  '"CRKSPH/SolidCRKSPHRZ.hh"',
                  '"CRKSPH/computeCRKSPHSumMassDensity.hh"',
                  '"CRKSPH/computeSolidCRKSPHSumMassDensity.hh"',
                  '"CRKSPH/detectSurface.hh"',
                  '"CRKSPH/centerOfMass.hh"',
                  '"CRKSPH/editMultimaterialSurfaceTopology.hh"',
                  '"CRKSPH/zerothOrderSurfaceCorrections.hh"',
                  '"Utilities/NodeCoupling.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
                  '"Neighbor/PairwiseField.hh"',
                  '"FileIO/FileIO.hh"',
                  '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeCRKSPHSumMassDensity(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                                W = "const TableKernel<%(Dimension)s>&",
                                position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                vol = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                                massDensity = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the CRKSPH mass density summation."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeSolidCRKSPHSumMassDensity(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                                     W = "const TableKernel<%(Dimension)s>&",
                                     position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                     mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                     H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                                     massDensity0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                     nodeCoupling = "const NodeCoupling&",
                                     massDensity = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the CRKSPH mass density summation."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def detectSurface(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                  m0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                  m1 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                  detectThreshold = "const double",
                  detectRange = "const double",
                  sweepAngle = "const double",
                  surfacePoint = "FieldList<%(Dimension)s, int>&"):
    "Detect surface particles leveraging the zeroth and first moments"
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def editMultimaterialSurfaceTopology(surfacePoint = "FieldList<%(Dimension)s, int>&",
                                     connectivityMap = "ConnectivityMap<%(Dimension)s>&"):
    """Look for any points that touch a surface (multi-material or void).
For such points:
  - Remove any non-surface multi-material overlap.
  - If not a surface point, flag this point as touching a surface point with
    surfacePoint=-1."""
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def zerothOrderSurfaceCorrections(A = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                  B = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                  C = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                                  gradA = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                  gradB = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                                  gradC = "FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                                  m0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                                  gradm0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                  surfacePoint = "const FieldList<%(Dimension)s, int>&"):
    """Look for any points that touch a surface (multi-material or void).
For such points:
  - enforce only zeroth order corrections."""
    return "void"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
CRKSPHBase%(ndim)id = PYB11TemplateClass(CRKSPHBase, template_parameters="%(Dimension)s")
CRKSPH%(ndim)id = PYB11TemplateClass(CRKSPH, template_parameters="%(Dimension)s")
SolidCRKSPH%(ndim)id = PYB11TemplateClass(SolidCRKSPH, template_parameters="%(Dimension)s")

@PYB11pycppname("centerOfMass")
def centerOfMass%(ndim)id(polyvol = "const %(Dimension)s::FacetedVolume&",
                          gradRhoi = "const %(Dimension)s::Vector&"):
    "Compute the center of mass of a FacetedVolume assuming a linear mass density field."
    return "%(Dimension)s::Vector"

computeCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
computeSolidCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeSolidCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
detectSurface%(ndim)id = PYB11TemplateFunction(detectSurface, template_parameters="%(Dimension)s")
editMultimaterialSurfaceTopology%(ndim)id = PYB11TemplateFunction(editMultimaterialSurfaceTopology, template_parameters="%(Dimension)s")
zerothOrderSurfaceCorrections%(ndim)id = PYB11TemplateFunction(zerothOrderSurfaceCorrections, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

#-------------------------------------------------------------------------------
# 2D
#-------------------------------------------------------------------------------
if 2 in dims:
    from CRKSPHRZ import *
    from SolidCRKSPHRZ import *
