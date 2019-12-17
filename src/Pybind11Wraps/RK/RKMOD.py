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
                  '"FileIO/FileIO.hh"',
                  '"Boundary/Boundary.hh"',
                  '"SPH/NodeCoupling.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

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
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
computeRKVolumes%(ndim)id = PYB11TemplateFunction(computeRKVolumes, template_parameters="%(Dimension)s")
computeRKSumVolume%(ndim)id = PYB11TemplateFunction(computeRKSumVolume, template_parameters="%(Dimension)s")
computeOccupancyVolume%(ndim)id = PYB11TemplateFunction(computeOccupancyVolume, template_parameters="%(Dimension)s")
computeVoronoiVolume%(ndim)id = PYB11TemplateFunction(computeVoronoiVolume, template_parameters="%(Dimension)s", pyname="computeVoronoiVolume")
computeHullVolumes%(ndim)id = PYB11TemplateFunction(computeHullVolumes, template_parameters="%(Dimension)s")
computeHVolumes%(ndim)id = PYB11TemplateFunction(computeHVolumes, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

    for num, correctionOrder in zip((0, 1, 2, 3, 4, 5, 6, 7),
                                    ("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder")):
        exec(''' 
RKCorrections%(ndim)sd%(num)s = PYB11TemplateClass(RKCorrections, template_parameters=("%(Dimension)s", "RKOrder::%(correctionOrder)s"))
RKUtilities%(ndim)sd%(num)s = PYB11TemplateClass(RKUtilities, template_parameters=("%(Dimension)s", "RKOrder::%(correctionOrder)s"))
''' % {"ndim"            : ndim,
       "Dimension"       : "Dim<" + str(ndim) + ">",
       "num"             : num,
       "correctionOrder" : correctionOrder})
