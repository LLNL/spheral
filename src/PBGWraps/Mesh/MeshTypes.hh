#ifndef __PBGWRAPS_MESHTYPES__
#define __PBGWRAPS_MESHTYPES__

#include "Geometry/Dimension.hh"
#include "Mesh/Mesh.hh"
#include "Mesh/Node.hh"
#include "Mesh/Edge.hh"
#include "Mesh/Face.hh"
#include "Mesh/Zone.hh"
#include "Mesh/computeGenerators.hh"
#include "Mesh/generateMesh.hh"
#include "Mesh/VoroPP.hh"
#include "Mesh/MeshWall.hh"
#include "Mesh/MeshConstructionUtilities.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {
namespace MeshSpace {

typedef Mesh<Dim<1> > LineMesh;
typedef Mesh<Dim<2> > PolygonalMesh;
typedef Mesh<Dim<3> > PolyhedralMesh;

typedef VoroPP<Dim<2> > VoroPP2d;
typedef VoroPP<Dim<3> > VoroPP3d;

typedef MeshWall<Dim<1> > MeshWall1d;
typedef MeshWall<Dim<2> > MeshWall2d;
typedef MeshWall<Dim<3> > MeshWall3d;

typedef PlanarMeshWall<Dim<1> > PlanarMeshWall1d;
typedef PlanarMeshWall<Dim<2> > PlanarMeshWall2d;
typedef PlanarMeshWall<Dim<3> > PlanarMeshWall3d;

typedef FacetedMeshWall<Dim<1> > FacetedMeshWall1d;
typedef FacetedMeshWall<Dim<2> > FacetedMeshWall2d;
typedef FacetedMeshWall<Dim<3> > FacetedMeshWall3d;

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to computeGenerators.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
computeGenerators(std::vector<NodeSpace::NodeList<Dimension>*>& nodeLists,
                  std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                  const typename Dimension::Vector& xmin,
                  const typename Dimension::Vector& xmax,
                  std::vector<typename Dimension::Vector>& positions,
                  std::vector<typename Dimension::SymTensor>& Hs,
                  std::vector<unsigned>& offsets) {
  computeGenerators<Dimension, 
                    typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator,
                    typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator>
    (nodeLists.begin(), nodeLists.end(), 
     boundaries.begin(), boundaries.end(),
     xmin, xmax, positions, Hs, offsets);
}

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to generateMesh
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
generateMesh(std::vector<NodeSpace::NodeList<Dimension>*>& nodeLists,
             std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
             const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const bool generateVoid,
             const bool generateParallelConnectivity,
             const bool removeBoundaryZones,
             const double voidThreshold,
             Mesh<Dimension>& mesh,
             NodeSpace::NodeList<Dimension>& voidNodes) {
  generateMesh(nodeLists.begin(), nodeLists.end(), 
               boundaries.begin(), boundaries.end(),
               xmin, xmax,
               generateVoid, generateParallelConnectivity, removeBoundaryZones, 
               voidThreshold,
               mesh, voidNodes);
}

//------------------------------------------------------------------------------
// Wrap the hashPosition method to return a Python tuple rather than a 
// boost::tuple.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PyObject*
hashPositionWrapper(const typename Dimension::Vector& position,
                    const typename Dimension::Vector& xmin,
                    const typename Dimension::Vector& xmax,
                    const typename Dimension::Vector& boxInv) {
  typename Mesh<Dimension>::Key result0 = hashPosition(position, xmin, xmax, boxInv);
  PyObject* result;
  result = Py_BuildValue("(KKK)", boost::get<0>(result0), boost::get<1>(result0), boost::get<2>(result0));
  return result;
}

}
}

#endif
