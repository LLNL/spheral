"""
Spheral Mesh module.

Provides useful utilities for defining and manipulating unstructured polygonal/polyhedral meshes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from Mesh import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Mesh/Mesh.hh"',
                  '"Mesh/Node.hh"',
                  '"Mesh/Edge.hh"',
                  '"Mesh/Face.hh"',
                  '"Mesh/Zone.hh"',
                  '"Mesh/computeGenerators.hh"',
                  '"Mesh/generateMesh.hh"',
                  '"Mesh/MeshConstructionUtilities.hh"',
                  '"FileIO/FileIO.hh"',
                  '"Boundary/Boundary.hh"',
                  '"NodeList/NodeList.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11implementation('''[](const std::vector<NodeList<%(Dimension)s>*>& nodeLists,
                           const std::vector<Boundary<%(Dimension)s>*>& boundaries,
                           const bool meshGhostNodes,
                           const typename %(Dimension)s::Vector& xmin,
                           const typename %(Dimension)s::Vector& xmax) {
                               std::vector<typename %(Dimension)s::Vector> positions;
                               std::vector<typename %(Dimension)s::SymTensor> Hs;
                               std::vector<unsigned> offsets;
                               computeGenerators<%(Dimension)s>(nodeLists.begin(), nodeLists.end(),
                                                                boundaries.begin(), boundaries.end(),
                                                                meshGhostNodes, xmin, xmax,
                                                                positions, Hs, offsets);
                               return py::make_tuple(positions, Hs, offsets);
                           }''')
def computeGenerators(nodeLists = "const std::vector<NodeList<%(Dimension)s>*>&",
                      boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                      meshGhostNodes = "const bool",
                      xmin = "const typename %(Dimension)s::Vector&",
                      xmax = "const typename %(Dimension)s::Vector&"):
    """This method takes a set of NodeLists, figures out what domains have to talk
to one another, and returns the flattened set of positions and Hs for the 
the generators this domain needs (including those from neighbor domains).
Return tuple contains (positions, Hs, offsets)."""
    return "py::tuple"

@PYB11template("Dimension")
@PYB11returnpolicy("take_ownership")
@PYB11implementation('''[](std::vector<NodeList<%(Dimension)s>*>& nodeLists,
                           std::vector<Boundary<%(Dimension)s>*>& boundaries,
                           const typename %(Dimension)s::Vector& xmin,
                           const typename %(Dimension)s::Vector& xmax,
                           const bool meshGhostNodes,
                           const bool generateVoid,
                           const bool generateParallelConnectivity,
                           const bool removeBoundaryZones,
                           const double voidThreshold,
                           Mesh<%(Dimension)s>& mesh,
                           NodeList<%(Dimension)s>& voidNodes) {
                               nodeLists.push_back(&voidNodes);
                               Spheral::generateMesh<%(Dimension)s>(nodeLists.begin(), nodeLists.end(),
                                                                    boundaries.begin(), boundaries.end(),
                                                                    xmin, xmax, meshGhostNodes, generateVoid,
                                                                    generateParallelConnectivity,
                                                                    removeBoundaryZones, voidThreshold,
                                                                    mesh, voidNodes);
                               nodeLists.pop_back();
                           }''')
def generateMesh(nodeLists = "const std::vector<NodeList<%(Dimension)s>*>&",
                 boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                 xmin = "const typename %(Dimension)s::Vector&",
                 xmax = "const typename %(Dimension)s::Vector&",
                 meshGhostNodes = "const bool",
                 generateVoid = "const bool",
                 generateParallelConnectivity = "const bool",
                 removeBoundaryZones = "const bool",
                 voidThreshold = "const double",
                 mesh = "Mesh<%(Dimension)s>&",
                 voidNodes = "NodeList<%(Dimension)s>&"):
    "Generate a mesh for the given set of NodeLists: returns tuple(mesh, voidNodes)"
    return None

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
%(prefix)sMesh = PYB11TemplateClass(Mesh, template_parameters="%(Dimension)s")

computeGenerators%(ndim)id = PYB11TemplateFunction(computeGenerators, template_parameters="%(Dimension)s")
generateMesh%(ndim)id = PYB11TemplateFunction(generateMesh, template_parameters="%(Dimension)s")

@PYB11pycppname("hashPosition")
@PYB11implementation("""[](const %(Vector)s& position,
                           const %(Vector)s& xmin,
                           const %(Vector)s& xmax,
                           const %(Vector)s& boxInv) -> py::tuple { 
    const auto result = hashPosition(position, xmin, xmax, boxInv);
    return py::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
  }""")
def hashPosition%(ndim)id(position = "const %(Vector)s&",
                          xmin = "const %(Vector)s&",
                          xmax = "const %(Vector)s&",
                          boxInv = "const %(Vector)s&"):
    "Hash a position for meshing -- complementary with quantizedPosition"
    return "py::tuple"

@PYB11pycppname("quantizedPosition")
@PYB11implementation("""[](const py::tuple hash, const %(Vector)s& xmin, const %(Vector)s& xmax) -> %(Vector)s { 
      if (hash.size() != 3) throw std::runtime_error("quantizedPosition ERROR: hash must be a 3-tuple");
      return quantizedPosition(std::make_tuple(hash[0].cast<uint64_t>(), hash[1].cast<uint64_t>(), hash[2].cast<uint64_t>()), xmin, xmax);
    }""")
def quantizedPosition%(ndim)id(hash = "py::tuple",
                               xmin = "const %(Vector)s&",
                               xmax = "const %(Vector)s&"):
    "Turn a position hash back into a cell center position -- complementary with hashPosition"
    return "%(Vector)s"
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "prefix"    : ("Line", "Polygonal", "Polyheral")[ndim-1],
       "Vector"    : "Dim<" + str(ndim) + ">::Vector"})
