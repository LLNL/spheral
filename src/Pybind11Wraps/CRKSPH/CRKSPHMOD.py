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
from RKCorrections import *
from RKUtilities import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"CRKSPH/CRKSPHUtilities.hh"',
                  '"CRKSPH/CRKSPHHydroBase.hh"',
                  '"CRKSPH/CRKSPHHydroBaseRZ.hh"',
                  '"CRKSPH/SolidCRKSPHHydroBase.hh"',
                  '"CRKSPH/SolidCRKSPHHydroBaseRZ.hh"',
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
                  '"CRKSPH/computeHVolumes.hh"',
                  '"CRKSPH/computeOccupancyVolume.hh"',
                  '"CRKSPH/gradientCRKSPH.hh"',
                  '"CRKSPH/interpolateCRKSPH.hh"',
                  '"CRKSPH/editMultimaterialSurfaceTopology.hh"',
                  '"CRKSPH/zerothOrderSurfaceCorrections.hh"',
                  '"CRKSPH/computeRKVolumes.hh"',
                  '"CRKSPH/RKCorrections.hh"',
                  '"CRKSPH/RKUtilities.hh"',
                  '"SPH/NodeCoupling.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
                  '<iterator>']

#-------------------------------------------------------------------------------
# Preamble
#-------------------------------------------------------------------------------
PYB11preamble += """
namespace {
template<typename Dimension>
struct AppendFieldLists: public boost::static_visitor<> {
    
    typedef std::vector<boost::variant<FieldList<Dimension, typename Dimension::Scalar>,
                                       FieldList<Dimension, typename Dimension::Vector>,
                                       FieldList<Dimension, typename Dimension::Tensor>,
                                       FieldList<Dimension, typename Dimension::SymTensor>,
                                       FieldList<Dimension, typename Dimension::ThirdRankTensor>>> FieldListArray;

    py::list& pylist;

    AppendFieldLists(py::list& pylist_):
        boost::static_visitor<>(),
        pylist(pylist_) {}

    template<typename FieldListType>
    inline
    void operator()(FieldListType& x) const {
        pylist.append(x);
    }
};
}
"""

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# enums
#-------------------------------------------------------------------------------
CRKOrder = PYB11enum(("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder"),
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

#-------------------------------------------------------------------------------
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
def computeCRKSPHSumVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           W = "const TableKernel<%(Dimension)s>&",
                           position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                           mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                           H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                           vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the CRKSPH mass density summation."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeOccupancyVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           W = "const TableKernel<%(Dimension)s>&",
                           position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                           H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                           vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the occupancy volume per point"
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
                             surfacePoint = "const FieldList<%(Dimension)s, int>&",
                             correctionOrder = "const CRKOrder",
                             A = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                             B = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             C = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                             gradA = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             gradB = "FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>&",
                             gradC = "FieldList<%(Dimension)s, typename %(Dimension)s::ThirdRankTensor>&"):
    "Function to compute CRK corrections based on the moments."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
@PYB11implementation("""[](const FieldList<%(Dimension)s, %(DataType)s>& fieldList,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& position,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& weight,
                           const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>& H,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& A,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& B,
                           const FieldList<%(Dimension)s, %(Dimension)s::Tensor>& C,
                           const ConnectivityMap<%(Dimension)s>& connectivityMap,
                           const CRKOrder correctionOrder,
                           const TableKernel<%(Dimension)s>& W,
                           const NodeCoupling& nodeCoupling) {
                               std::vector<boost::variant<FieldList<%(Dimension)s, %(Dimension)s::Scalar>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::Vector>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::Tensor>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::SymTensor>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>>> fieldLists;
                               fieldLists.emplace_back(fieldList);
                               auto flvec = interpolateCRKSPH(fieldLists,
                                                              position,
                                                              weight,
                                                              H,
                                                              A,
                                                              B,
                                                              C,
                                                              connectivityMap,
                                                              correctionOrder,
                                                              W,
                                                              nodeCoupling);
                               CHECK(flvec.size() == 1);
                               FieldList<%(Dimension)s, %(DataType)s> result(boost::get<FieldList<%(Dimension)s, %(DataType)s>>(flvec[0]));
                               result.copyFields();
                               return result;
                           }""")
@PYB11cppname("interpolateCRKSPH")
def interpolateCRKSPH1(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
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

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11implementation("""[](const py::list& pyFieldLists,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& position,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& weight,
                           const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>& H,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& A,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& B,
                           const FieldList<%(Dimension)s, %(Dimension)s::Tensor>& C,
                           const ConnectivityMap<%(Dimension)s>& connectivityMap,
                           const CRKOrder correctionOrder,
                           const TableKernel<%(Dimension)s>& W,
                           const NodeCoupling& nodeCoupling) {
                               std::vector<boost::variant<FieldList<%(Dimension)s, %(Dimension)s::Scalar>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::Vector>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::Tensor>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::SymTensor>,
                                                          FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>>> fieldLists;
                               for (auto& pyfl: pyFieldLists) {
                                   try {
                                           auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Scalar>&>();
                                           fieldLists.emplace_back(fl);
                                   } catch (py::cast_error& e) {
                                       try {
                                               auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Vector>&>();
                                               fieldLists.emplace_back(fl);
                                       } catch (py::cast_error& e) {
                                           try {
                                                   auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Tensor>&>();
                                                   fieldLists.emplace_back(fl);
                                           } catch (py::cast_error& e) {
                                               try {
                                                       auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&>();
                                                       fieldLists.emplace_back(fl);
                                               } catch (py::cast_error& e) {
                                                   try {
                                                           auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>&>();
                                                           fieldLists.emplace_back(fl);
                                                   } catch (py::cast_error& e) {
                                                       throw py::type_error("interpolateCRKSPH expects a list of only (Scalar,Vector,Tensor,SymTensor,ThirdRankTensor)FieldList.");
                                                   }
                                               }
                                           }
                                       }
                                   }
                               }
                               auto cppresult =  interpolateCRKSPH(fieldLists,
                                                                   position,
                                                                   weight,
                                                                   H,
                                                                   A,
                                                                   B,
                                                                   C,
                                                                   connectivityMap,
                                                                   correctionOrder,
                                                                   W,
                                                                   nodeCoupling);
                               py::list result;
                               for (auto fl: cppresult) boost::apply_visitor(AppendFieldLists<%(Dimension)s>(result), fl);
                               return result;
                           }""")
def interpolateCRKSPH(fieldLists = "py::list&",
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
    return "py::list"

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeVoronoiVolume(position = "const FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
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
                         cellFaceFlags = "FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&"):
    "Compute the volume per point based on the Voronoi tessellation-like algorithm."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeHullVolumes(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                       kernelExtent = "const typename %(Dimension)s::Scalar",
                       position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                       H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                       volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the volume per point based on convex hulls."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeHVolumes(nPerh = "const typename %(Dimension)s::Scalar",
                    H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                    volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute a volume per point based on the local H tensor."
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
@PYB11template("Dimension")
def computeRKVolumes(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                     W = "const TableKernel<%(Dimension)s>&",
                     position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                     mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                     massDensity = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                     H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                     damage = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                     boundaryConditions = "const std::vector<Boundary<%(Dimension)s>*>&",
                     volumeType = "const CRKVolumeType",
                     surfacePoint = "FieldList<%(Dimension)s, int>&",
                     deltaCentroid = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                     etaVoidPoints = "const FieldList<%(Dimension)s, std::vector<typename %(Dimension)s::Vector>>&",
                     cells = "FieldList<%(Dimension)s, typename %(Dimension)s::FacetedVolume>&",
                     cellFaceFlags = "FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&",
                     volume = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute RK volumes"
    return "void"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
CRKSPHHydroBase%(ndim)id = PYB11TemplateClass(CRKSPHHydroBase, template_parameters="%(Dimension)s")
SolidCRKSPHHydroBase%(ndim)id = PYB11TemplateClass(SolidCRKSPHHydroBase, template_parameters="%(Dimension)s")

@PYB11pycppname("centerOfMass")
def centerOfMass%(ndim)id(polyvol = "const %(Dimension)s::FacetedVolume&",
                          gradRhoi = "const %(Dimension)s::Vector&"):
    "Compute the center of mass of a FacetedVolume assuming a linear mass density field."
    return "%(Dimension)s::Vector"

CRKSPHKernel%(ndim)id = PYB11TemplateFunction(CRKSPHKernel, template_parameters="%(Dimension)s")
CRKSPHKernelAndGradient%(ndim)id = PYB11TemplateFunction(CRKSPHKernelAndGradient, template_parameters="%(Dimension)s")
computeCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
computeCRKSPHSumVolume%(ndim)id = PYB11TemplateFunction(computeCRKSPHSumVolume, template_parameters="%(Dimension)s")
computeVoronoiVolume%(ndim)id = PYB11TemplateFunction(computeVoronoiVolume, template_parameters="%(Dimension)s", pyname="computeVoronoiVolume")
computeOccupancyVolume%(ndim)id = PYB11TemplateFunction(computeOccupancyVolume, template_parameters="%(Dimension)s")
computeSolidCRKSPHSumMassDensity%(ndim)id = PYB11TemplateFunction(computeSolidCRKSPHSumMassDensity, template_parameters="%(Dimension)s")
computeCRKSPHMoments%(ndim)id = PYB11TemplateFunction(computeCRKSPHMoments, template_parameters="%(Dimension)s")
detectSurface%(ndim)id = PYB11TemplateFunction(detectSurface, template_parameters="%(Dimension)s")
computeCRKSPHCorrections%(ndim)id = PYB11TemplateFunction(computeCRKSPHCorrections, template_parameters="%(Dimension)s")
computeHullVolumes%(ndim)id = PYB11TemplateFunction(computeHullVolumes, template_parameters="%(Dimension)s")
computeHVolumes%(ndim)id = PYB11TemplateFunction(computeHVolumes, template_parameters="%(Dimension)s")
editMultimaterialSurfaceTopology%(ndim)id = PYB11TemplateFunction(editMultimaterialSurfaceTopology, template_parameters="%(Dimension)s")
zerothOrderSurfaceCorrections%(ndim)id = PYB11TemplateFunction(zerothOrderSurfaceCorrections, template_parameters="%(Dimension)s")
computeRKVolumes%(ndim)id = PYB11TemplateFunction(computeRKVolumes, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

    # CRKSPH interpolation
    exec('''interpolateCRKSPH%(ndim)id = PYB11TemplateFunction(interpolateCRKSPH, template_parameters="Dim<%(ndim)i>", pyname="interpolateCRKSPH")''' % {"ndim" : ndim})
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim,
                    "Dim<%i>::Tensor" % ndim,
                    "Dim<%i>::SymTensor" % ndim,
                    "Dim<%i>::ThirdRankTensor" % ndim):
        exec('''
interpolateCRKSPH%(label)s = PYB11TemplateFunction(interpolateCRKSPH1, template_parameters=("%(Dimension)s", "%(element)s"), pyname="interpolateCRKSPH")
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

    # New RK classes
    for num, correctionOrder in zip((0, 1, 2, 3, 4, 5, 6, 7),
                                    ("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder")):
        exec(''' 
RKCorrections%(ndim)sd%(num)s = PYB11TemplateClass(RKCorrections, template_parameters=("%(Dimension)s", "CRKOrder::%(correctionOrder)s"))
RKUtilities%(ndim)sd%(num)s = PYB11TemplateClass(RKUtilities, template_parameters=("%(Dimension)s", "CRKOrder::%(correctionOrder)s"))
''' % {"ndim"            : ndim,
       "Dimension"       : "Dim<" + str(ndim) + ">",
       "num"             : num,
       "correctionOrder" : correctionOrder})

#-------------------------------------------------------------------------------
# 2D
#-------------------------------------------------------------------------------
if 2 in dims:
    from CRKSPHHydroBaseRZ import *
    from SolidCRKSPHHydroBaseRZ import *
