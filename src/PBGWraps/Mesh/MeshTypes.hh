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
#include "Mesh/MeshConstructionUtilities.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {

typedef Mesh<Dim<1> > LineMesh;
typedef Mesh<Dim<2> > PolygonalMesh;
typedef Mesh<Dim<3> > PolyhedralMesh;

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to computeGenerators.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
computeGenerators(std::vector<NodeList<Dimension>*>& nodeLists,
                  std::vector<Boundary<Dimension>*>& boundaries,
                  const bool meshGhostNodes,
                  const typename Dimension::Vector& xmin,
                  const typename Dimension::Vector& xmax,
                  std::vector<typename Dimension::Vector>& positions,
                  std::vector<typename Dimension::SymTensor>& Hs,
                  std::vector<unsigned>& offsets) {
  computeGenerators<Dimension, 
                    typename std::vector<NodeList<Dimension>*>::const_iterator,
                    typename std::vector<Boundary<Dimension>*>::const_iterator>
    (nodeLists.begin(), nodeLists.end(), 
     boundaries.begin(), boundaries.end(),
     meshGhostNodes, xmin, xmax, 
     positions, Hs, offsets);
}

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to generateMesh
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
generateMesh(std::vector<NodeList<Dimension>*>& nodeLists,
             std::vector<Boundary<Dimension>*>& boundaries,
             const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const bool meshGhostNodes,
             const bool generateVoid,
             const bool generateParallelConnectivity,
             const bool removeBoundaryZones,
             const double voidThreshold,
             Mesh<Dimension>& mesh,
             NodeList<Dimension>& voidNodes) {
  generateMesh(nodeLists.begin(), nodeLists.end(), 
               boundaries.begin(), boundaries.end(),
               xmin, xmax,
               meshGhostNodes, generateVoid, generateParallelConnectivity, removeBoundaryZones, 
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

//------------------------------------------------------------------------------
// Similarly wrap the inverse operation, quantizedPosition.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
quantizedPositionWrapper(PyObject* hashTup,
                         const typename Dimension::Vector& xmin,
                         const typename Dimension::Vector& xmax) {
  VERIFY2(PyTuple_Check(hashTup) == true,
          "quantizedPosition ERROR: first argument must a tuple.");
  VERIFY2(PyTuple_Size(hashTup) == 3,
          "quantizedPosition ERROR: tuple must of length 3.");
  unsigned long long ix, iy, iz;
  VERIFY2(PyArg_ParseTuple(hashTup, (char*)"KKK", &ix, &iy, &iz),
          "quantizedPosition ERROR: unable to parse input tuple.");
  return quantizedPosition(boost::tuple<uint64_t, uint64_t, uint64_t>(ix, iy, iz), xmin, xmax);
}

//------------------------------------------------------------------------------
// Wrap the offset lookup up method to return a tuple -- nicer python interface.
//------------------------------------------------------------------------------
template<typename MeshType>
inline
PyObject*
lookupNodeListID(MeshType* self,
                 unsigned zoneID) {
  unsigned nodeListi, i;
  self->lookupNodeListID(zoneID, nodeListi, i);
  PyObject* result;
  result = Py_BuildValue("(II)", nodeListi, i);
  return result;
}

}

#endif
