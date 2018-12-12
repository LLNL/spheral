"""
Spheral CXXTests module.

This module provide thin front-end wrappers for the C++ unit test in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

PYB11includes = ['"CXXTests/testNodeIterators.hh"',
                 '"CXXTests/test_r3d_utils.hh"',
                 '"Geometry/Dimension.hh"',
                 '"DataBase/DataBase.hh"']

PYB11namespaces = ["Spheral"]

# #-------------------------------------------------------------------------------
# # Node iterator tests
# #-------------------------------------------------------------------------------
# @PYB11template("Dimension")
# def testGlobalAllNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global AllNodeIterators."
#     return "std::string"

# @PYB11template("Dimension")
# def testGlobalInternalNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global InternalNodeIterators."
#     return "std::string"

# @PYB11template("Dimension")
# def testGlobalGhostNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global GhostNodeIterators."
#     return "std::string"

# @PYB11template("Dimension")
# def testGlobalMasterNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global MasterNodeIterators."
#     return "std::string"

# @PYB11template("Dimension")
# def testGlobalCoarseNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global CoarseNodeIterators."
#     return "std::string"

# @PYB11template("Dimension")
# def testGlobalRefineNodeIterators(dataBase = "const DataBase<Dimension>&"):
#     "Test global RefineNodeIterators."
#     return "std::string"

# for ndim in dims:
#     exec('''
# testGlobalAllNodeIterators%(ndim)id      = PYB11TemplateFunction(testGlobalAllNodeIterators,      template_parameters="Dim<%(ndim)i>", pyname="testGlobalAllNodeIterators")
# testGlobalInternalNodeIterators%(ndim)id = PYB11TemplateFunction(testGlobalInternalNodeIterators, template_parameters="Dim<%(ndim)i>", pyname="testGlobalInternalNodeIterators")
# testGlobalGhostNodeIterators%(ndim)id    = PYB11TemplateFunction(testGlobalGhostNodeIterators,    template_parameters="Dim<%(ndim)i>", pyname="testGlobalGhostNodeIterators")
# testGlobalMasterNodeIterators%(ndim)id   = PYB11TemplateFunction(testGlobalMasterNodeIterators,   template_parameters="Dim<%(ndim)i>", pyname="testGlobalMasterNodeIterators")
# testGlobalCoarseNodeIterators%(ndim)id   = PYB11TemplateFunction(testGlobalCoarseNodeIterators,   template_parameters="Dim<%(ndim)i>", pyname="testGlobalCoarseNodeIterators")
# testGlobalRefineNodeIterators%(ndim)id   = PYB11TemplateFunction(testGlobalRefineNodeIterators,   template_parameters="Dim<%(ndim)i>", pyname="testGlobalRefineNodeIterators")
# ''' % {"ndim" : ndim})

#-------------------------------------------------------------------------------
# R3D tests
#-------------------------------------------------------------------------------
def test_polygon_to_r2d_poly():
    "Test converting r2d_poly <- polygon."
    return "std::string"

def test_r2d_poly_to_polygon():
    "Test converting r2d_poly -> polygon."
    return "std::string"

def test_polyhedron_to_r3d_poly():
    "Test converting r3d_poly <- polyhedron."
    return "std::string"

def test_r3d_poly_to_polyhedron():
    "Test converting r3d_poly -> polyhedron."
    return "std::string"

def test_clip_polygon():
    "Test clipping."
    return "std::string"

def test_orphan_polygon():
    "Test clipping."
    return "std::string"

def test_clip_polyhedron():
    "Test clipping."
    return "std::string"

