"""
Spheral RK module.

Provides reproducing kernel infrastructure.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from RKCorrections import *
from RKUtilities import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"RK/RKCorrections.hh"',
                  '"RK/RKUtilities.hh"',
                  '"RK/computeRKVolumes.hh"',
                  '"RK/computeVoronoiVolume.hh"',
                  '"RK/computeOccupancyVolume.hh"',
                  '"RK/computeRKSumVolume.hh"',
                  '"RK/computeHullVolumes.hh"',
                  '"RK/computeHVolumes.hh"',
                  '"RK/computeOccupancyVolume.hh"',
                  '"RK/interpolateRK.hh"',
                  '"FileIO/FileIO.hh"',
                  '"Boundary/Boundary.hh"',
                  '"SPH/NodeCoupling.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

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
# enums
#-------------------------------------------------------------------------------
RKOrder = PYB11enum(("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder"),
                    export_values = True,
                    doc = "Selection of RK correction orders")
RKVolumeType = PYB11enum(("RKMassOverDensity", "RKSumVolume", "RKVoronoiVolume", "RKHullVolume", "HVolume"),
                         export_values = True,
                         doc = "Options for RK volume algorithms")

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def RKKernel(order = "const RKOrder",
             W = "const TableKernel<%(Dimension)s>&",
             x = "const typename %(Dimension)s::Vector&",
             H = "const typename %(Dimension)s::SymTensor&",
             corrections = "const std::vector<double>&"):
    "Compute the corrected kernel value."
    return "typename %(Dimension)s::Scalar"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def RKGradient(order = "const RKOrder",
               W = "const TableKernel<%(Dimension)s>&",
               x = "const typename %(Dimension)s::Vector&",
               H = "const typename %(Dimension)s::SymTensor&",
               corrections = "const std::vector<double>&"):
    "Compute the corrected kernel gradient."
    return "typename %(Dimension)s::Vector"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11implementation("""[](const RKOrder order,
                           const TableKernel<%(Dimension)s>& W,
                           const typename %(Dimension)s::Vector& x,
                           const typename %(Dimension)s::SymTensor& H,
                           const std::vector<double>& corrections) {
                               typename %(Dimension)s::Scalar WRK, gradWSPH;
                               typename %(Dimension)s::Vector gradWRK;
                               RKKernelAndGradient(WRK, gradWSPH, gradWRK,                            
                                                   order,
                                                   W,
                                                   x,
                                                   H,
                                                   corrections);
                               return py::make_tuple(WRK, gradWSPH, gradWRK);
                           }""")
def RKKernelAndGradient(order = "const RKOrder",
                        W = "const TableKernel<%(Dimension)s>&",
                        x = "const typename %(Dimension)s::Vector&",
                        H = "const typename %(Dimension)s::SymTensor&",
                        corrections = "const std::vector<double>&"):
    "Compute the corrected kernel, SPH gradient magnitude, and corrected gradient."
    return "py::tuple"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeRKCorrections(order = "const RKOrder",
                         connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                         W = "const TableKernel<%(Dimension)s>&",
                         volume = "const FieldList<%(Dimension)s, %(Scalar)s>&",
                         position = "const FieldList<%(Dimension)s, %(Vector)s>&",
                         H = "const FieldList<%(Dimension)s, %(SymTensor)s>&",
                         needHessian = "const bool",
                         zerothCorrections = "FieldList<%(Dimension)s, std::vector<double>>&",
                         corrections = "FieldList<%(Dimension)s, std::vector<double>>&"):
    "Compute the RK corrections for the requested order"
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
                     volumeType = "const RKVolumeType",
                     surfacePoint = "FieldList<%(Dimension)s, int>&",
                     deltaCentroid = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                     etaVoidPoints = "const FieldList<%(Dimension)s, std::vector<typename %(Dimension)s::Vector>>&",
                     cells = "FieldList<%(Dimension)s, typename %(Dimension)s::FacetedVolume>&",
                     cellFaceFlags = "FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&",
                     volume = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute RK volumes"
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeRKSumVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                       W = "const TableKernel<%(Dimension)s>&",
                       position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                       mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                       H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                       vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the RK mass density summation."
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
@PYB11implementation("""[](const py::list& pyFieldLists,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& position,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& weight,
                           const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>& H,
                           const ConnectivityMap<%(Dimension)s>& connectivityMap,
                           const TableKernel<%(Dimension)s>& W,
                           const RKOrder correctionOrder,
                           const FieldList<%(Dimension)s, std::vector<double>>& corrections,
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
                               auto cppresult =  interpolateRK(fieldLists,
                                                               position,
                                                               weight,
                                                               H,
                                                               connectivityMap,
                                                               W,
                                                               correctionOrder,
                                                               corrections,
                                                               nodeCoupling);
                               py::list result;
                               for (auto fl: cppresult) boost::apply_visitor(AppendFieldLists<%(Dimension)s>(result), fl);
                               return result;
                           }""")
def interpolateRK(fieldLists = "py::list&",
                  position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                  H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                  connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                  W = "const TableKernel<%(Dimension)s>&",
                  correctionOrder = "const RKOrder",
                  corrections = "const FieldList<%(Dimension)s, std::vector<double>>&",
                  nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK interpolation at each point."
    return "py::list"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
RKKernel%(ndim)id = PYB11TemplateFunction(RKKernel, template_parameters="%(Dimension)s")
RKGradient%(ndim)id = PYB11TemplateFunction(RKGradient, template_parameters="%(Dimension)s")
RKKernelAndGradient%(ndim)id = PYB11TemplateFunction(RKKernelAndGradient, template_parameters="%(Dimension)s")
computeRKCorrections%(ndim)id = PYB11TemplateFunction(computeRKCorrections, template_parameters="%(Dimension)s")
computeRKVolumes%(ndim)id = PYB11TemplateFunction(computeRKVolumes, template_parameters="%(Dimension)s")
computeRKSumVolume%(ndim)id = PYB11TemplateFunction(computeRKSumVolume, template_parameters="%(Dimension)s")
computeOccupancyVolume%(ndim)id = PYB11TemplateFunction(computeOccupancyVolume, template_parameters="%(Dimension)s")
computeVoronoiVolume%(ndim)id = PYB11TemplateFunction(computeVoronoiVolume, template_parameters="%(Dimension)s", pyname="computeVoronoiVolume")
computeHullVolumes%(ndim)id = PYB11TemplateFunction(computeHullVolumes, template_parameters="%(Dimension)s")
computeHVolumes%(ndim)id = PYB11TemplateFunction(computeHVolumes, template_parameters="%(Dimension)s")
interpolateRK%(ndim)id = PYB11TemplateFunction(interpolateRK, template_parameters="Dim<%(ndim)i>", pyname="interpolateRK")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "Scalar"    : "Dim<" + str(ndim) + ">::Scalar",
       "Vector"    : "Dim<" + str(ndim) + ">::Vector",
       "Tensor"    : "Dim<" + str(ndim) + ">::Tensor",
       "SymTensor" : "Dim<" + str(ndim) + ">::SymTensor"})

    # RKCorrections and RKUtilities
    for num, correctionOrder in enumerate(("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder")):
        exec(''' 
RKCorrections%(ndim)sd%(num)s = PYB11TemplateClass(RKCorrections, template_parameters=("%(Dimension)s", "RKOrder::%(correctionOrder)s"))
RKUtilities%(ndim)sd%(num)s = PYB11TemplateClass(RKUtilities, template_parameters=("%(Dimension)s", "RKOrder::%(correctionOrder)s"))
''' % {"ndim"            : ndim,
       "Dimension"       : "Dim<" + str(ndim) + ">",
       "num"             : num,
       "correctionOrder" : correctionOrder})

    # RK interpolation
