"""
Spheral CRKSPH module.

Provides implementations of CRKSPH (fluid and solid).
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from CRKSPHHydroBase import *
from SolidCRKSPHHydroBase import *
from CRKSPHVariant import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"CRKSPH/CRKSPHUtilities.hh"',
                  '"CRKSPH/CRKSPHHydroBase.hh"',
                  '"CRKSPH/CRKSPHHydroBaseRZ.hh"',
                  '"CRKSPH/SolidCRKSPHHydroBase.hh"',
                  '"CRKSPH/SolidCRKSPHHydroBaseRZ.hh"',
                  '"CRKSPH/CRKSPHVariant.hh"',
                  '"CRKSPH/computeVoronoiVolume.hh"',
                  '"CRKSPH/computeOccupancyVolume.hh"',
                  '"CRKSPH/computeCRKSPHSumVolume.hh"',
                  '"CRKSPH/computeCRKSPHSumMassDensity.hh"',
                  '"CRKSPH/computeSolidCRKSPHSumMassDensity.hh"',
                  '"CRKSPH/computeCRKSPHMoments.hh"',
                  '"CRKSPH/detectSurface.hh"',
                  '"CRKSPH/computeCRKSPHCorrections.hh"',
                  '"CRKSPH/centerOfMass.hh"',
                  '"CRKSPH/computeHullVolumes.hh"',
                  '"CRKSPH/computeNeighborHull.hh"',
                  '"CRKSPH/computeHVolumes.hh"',
                  '"CRKSPH/computeOccupancyVolume.hh"',
                  '"CRKSPH/gradientCRKSPH.hh"',
                  '"CRKSPH/interpolateCRKSPH.hh"',
                  '"SPH/NodeCoupling.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
                  '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# enums
#-------------------------------------------------------------------------------
CRKOrder = PYB11enum(("ZerothOrder", "LinearOrder", "QuadraticOrder"),
                     export_values = True,
                     doc = "Selection of CRK correction orders")
CRKVolumeType = PYB11enum(("CRKMassOverDensity", "CRKSumVolume", "CRKVoronoiVolume", "CRKHullVolume", "HVolume"),
                          export_values = True,
                          doc = "Options for CRK mass density algorithms")

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def CRKSPHKernel(W = "const TableKernel<%(Dimension)s>&",
                 correctionOrder = "const CRKOrder",
                 rij = "const typename %(Dimension)s::Vector&",
                 etaj = "const typename %(Dimension)s::Vector&",
                 Hdetj = "const typename %(Dimension)s::Scalar",
                 Ai = "const typename %(Dimension)s::Scalar",
                 Bi = "const typename %(Dimension)s::Vector&",
                 Ci = "const typename %(Dimension)s::Tensor&"):
    "Compute the corrected kernel value."
    return "typename %(Dimension)s::Scalar"

@PYB11template("Dimension")
@PYB11implementation("""[](const TableKernel<%(Dimension)s>& W,
                           const CRKOrder correctionOrder,
                           const typename %(Dimension)s::Vector& rij,
                           const typename %(Dimension)s::Vector& etaj,
                           const typename %(Dimension)s::SymTensor& Hj,
                           const typename %(Dimension)s::Scalar Hdetj,
                           const typename %(Dimension)s::Scalar Ai,
                           const typename %(Dimension)s::Vector& Bi,
                           const typename %(Dimension)s::Tensor& Ci,
                           const typename %(Dimension)s::Vector& gradAi,
                           const typename %(Dimension)s::Tensor& gradBi,
                           const typename %(Dimension)s::ThirdRankTensor& gradCi) {
                               typename %(Dimension)s::Scalar WCRKSPH, gradWSPH;
                               typename %(Dimension)s::Vector gradWCRKSPH;
                               CRKSPHKernelAndGradient(WCRKSPH, gradWSPH, gradWCRKSPH,                            
                                                       W,
                                                       correctionOrder,
                                                       rij,
                                                       etaj,
                                                       Hj,
                                                       Hdetj,
                                                       Ai,
                                                       Bi,
                                                       Ci,
                                                       gradAi,
                                                       gradBi,
                                                       gradCi);
                               return py::make_tuple(WCRKSPH, gradWSPH, gradWCRKSPH);
                           }""")
def CRKSPHKernelAndGradient(
        W = "const TableKernel<%(Dimension)s>&",
        correctionOrder = "const CRKOrder",
        rij = "const typename %(Dimension)s::Vector&",
        etaj = "const typename %(Dimension)s::Vector&",
        Hj = "const typename %(Dimension)s::SymTensor&",
        Hdetj = "const typename %(Dimension)s::Scalar",
        Ai = "const typename %(Dimension)s::Scalar",
        Bi = "const typename %(Dimension)s::Vector&",
        Ci = "const typename %(Dimension)s::Tensor&",
        gradAi = "const typename %(Dimension)s::Vector&",
        gradBi = "const typename %(Dimension)s::Tensor&",
        gradCi = "const typename %(Dimension)s::ThirdRankTensor&"):
    "Compute the corrected kernel value, uncorrected and corrected gradients."
    return "py::tuple"

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

@PYB11template("Dimension")
def computeCRKSPHSumVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           W = "const TableKernel<%(Dimension)s>&",
                           position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                           mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                           H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                           vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the CRKSPH mass density summation."
    return "void"

@PYB11template("Dimension")
def computeOccupancyVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           W = "const TableKernel<%(Dimension)s>&",
                           position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                           H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                           vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the occupancy volume per point"
    return "void"

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

@PYB11template("Dimension")
def computeCRKSPHMoments(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                         W = "const TableKernel<%(Dimension)s>&",
                         weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                         position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                         H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                         correctionOrder = "const CRKOrder",
                         nodeCoupling = "const NodeCoupling&",
                         m0 = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                         m1 = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                         m2 = "FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                         m3 = "FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                         m4 = "FieldList<%(Dimension)s, typename %(Dimension)s::FourthRankTensor>&",
                         gradm0 = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                         gradm1 = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                         gradm2 = "FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                         gradm3 = "FieldList<%(Dimension)s, typename %(Dimension)s::FourthRankTensor>&",
                         gradm4 = "FieldList<%(Dimension)s, typename %(Dimension)s::FifthRankTensor>&"):
    "Compute the moments necessary for CRKSPH corrections."
    return "void"

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

@PYB11template("Dimension")
def computeCRKSPHCorrections(m0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                             m1 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             m2 = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                             m3 = "const FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                             m4 = "const FieldList<%(Dimension)s, typename %(Dimension)s::FourthRankTensor>&",
                             gradm0 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             gradm1 = "const FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                             gradm2 = "const FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                             gradm3 = "const FieldList<%(Dimension)s, typename %(Dimension)s::FourthRankTensor>&",
                             gradm4 = "const FieldList<%(Dimension)s, typename %(Dimension)s::FifthRankTensor>&",
                             H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                             correctionOrder = "const CRKOrder",
                             A = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                             B = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             C = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                             gradA = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             gradB = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                             gradC = "FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&"):
    "Function to compute CRK corrections based on the moments."
    return "void"

@PYB11template("Dimension", "DataType")
def interpolateCRKSPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                      position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                      weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                      H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                      A = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                      B = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                      C = "const FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                      connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                      correctionOrder = "const CRKOrder",
                      W = "const TableKernel<%(Dimension)s>&",
                      nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the CRKSPH interpolation at each point."
    return "FieldList<%(Dimension)s, %(DataType)s>"

@PYB11template("Dimension", "DataType")
def gradientCRKSPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                   position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                   weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                   H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                   A = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                   B = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                   C = "const FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                   gradA = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                   gradB = "const FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                   gradC = "const FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&",
                   connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                   correctionOrder = "const CRKOrder",
                   W = "const TableKernel<%(Dimension)s>&",
                   nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the CRKSPH gradient."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>"

@PYB11template("Dimension")
def computeHullVolumes(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                       kernelExtent = "const typename %(Dimension)s::Scalar",
                       position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                       H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                       volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the volume per point based on convex hulls."
    return "void"

@PYB11template("Dimension")
def computeHVolumes(nPerh = "const typename %(Dimension)s::Scalar",
                    H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                    volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute a volume per point based on the local H tensor."
    return "void"

@PYB11template("Dimension")
def computeNeighborHull(fullConnectivity = "const std::vector<std::vector<int> >&",
                        etaCutoff = "const typename %(Dimension)s::Scalar",
                        ri = "const typename %(Dimension)s::Vector&",
                        Hi = "const typename %(Dimension)s::SymTensor&",
                        position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&"):
    "Compute the hull for a given points neighbor set."
    return "typename %(Dimension)s::FacetedVolume"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
CRKSPHHydroBase%(ndim)id = PYB11TemplateClass(CRKSPHHydroBase, template_parameters="%(Dimension)s")
SolidCRKSPHHydroBase%(ndim)id = PYB11TemplateClass(SolidCRKSPHHydroBase, template_parameters="%(Dimension)s")
CRKSPHVariant%(ndim)id = PYB11TemplateClass(CRKSPHVariant, template_parameters="%(Dimension)s")

@PYB11pycppname("computeVoronoiVolume")
def computeVoronoiVolume%(ndim)id(position = "const FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
                                  H = "const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&",
                                  connectivityMap = "const ConnectivityMap<%(Dimension)s >&",
                                  damage = "const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&",
                                  facetedBoundaries = "const std::vector<%(Dimension)s::FacetedVolume>&",
                                  holes = "const std::vector<std::vector<%(Dimension)s::FacetedVolume> >&",
                                  boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                                  weight = "const FieldList<%(Dimension)s, %(Dimension)s::Scalar>&",
                                  surfacePoint = "FieldList<%(Dimension)s, int>&",
                                  vol = "FieldList<%(Dimension)s, %(Dimension)s::Scalar>&",
                                  deltaMedian = "FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
                                  etaVoidPoints = "FieldList<%(Dimension)s, std::vector<%(Dimension)s::Vector>>&",
                                  cells = "FieldList<%(Dimension)s, %(Dimension)s::FacetedVolume>&",
                                  cellFaceFlags = "FieldList<%(Dimension)s, std::vector<int>>&"):
    "Compute the volume per point based on the Voronoi tessellation-like algorithm."
    return "void"

@PYB11pycppname("centerOfMass")
def centerOfMass%(ndim)id(polyvol = "const %(Dimension)s::FacetedVolume&",
                          gradRhoi = "const %(Dimension)s::Vector&"):
    "Compute the center of mass of a FacetedVolume assuming a linear mass density field."
    return "%(Dimension)s::Vector"

CRKSPHKernel%(ndim)id = PYB11TemplateFunction(CRKSPHKernel, template_parameters="%(Dimension)s")
CRKSPHKernelAndGradient%(ndim)id = PYB11TemplateFunction(CRKSPHKernelAndGradient, template_parameters="%(Dimension)s")
computeCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
computeCRKSPHSumVolume%(ndim)id = PYB11TemplateFunction(computeCRKSPHSumVolume, template_parameters="%(Dimension)s")
computeOccupancyVolume%(ndim)id = PYB11TemplateFunction(computeOccupancyVolume, template_parameters="%(Dimension)s")
computeSolidCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeSolidCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
computeCRKSPHMoments%(ndim)id = PYB11TemplateFunction(computeCRKSPHMoments, template_parameters="%(Dimension)s")
detectSurface%(ndim)id = PYB11TemplateFunction(detectSurface, template_parameters="%(Dimension)s")
computeCRKSPHCorrections%(ndim)id = PYB11TemplateFunction(computeCRKSPHCorrections, template_parameters="%(Dimension)s")
computeHullVolumes%(ndim)id = PYB11TemplateFunction(computeHullVolumes, template_parameters="%(Dimension)s")
computeHVolumes%(ndim)id = PYB11TemplateFunction(computeHVolumes, template_parameters="%(Dimension)s")
computeNeighborHull%(ndim)id = PYB11TemplateFunction(computeNeighborHull, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

    # CRKSPH interpolation
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim,
                    "Dim<%i>::Tensor" % ndim,
                    "Dim<%i>::SymTensor" % ndim,
                    "Dim<%i>::ThirdRankTensor" % ndim):
        exec('''
interpolateCRKSPH%(label)s = PYB11TemplateFunction(interpolateCRKSPH, template_parameters=("%(Dimension)s", "%(element)s"), pyname="interpolateCRKSPH")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

    # CRKSPH gradient
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim):
        exec('''
gradientCRKSPH%(label)s = PYB11TemplateFunction(gradientCRKSPH, template_parameters=("%(Dimension)s", "%(element)s"), pyname="gradientCRKSPH")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

#-------------------------------------------------------------------------------
# 2D
#-------------------------------------------------------------------------------
if 2 in dims:
    from CRKSPHHydroBaseRZ import *
    from SolidCRKSPHHydroBaseRZ import *
