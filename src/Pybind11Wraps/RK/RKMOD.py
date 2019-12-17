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
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
computeRKVolumes%(ndim)id = PYB11TemplateFunction(computeRKVolumes, template_parameters="%(Dimension)s")
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
