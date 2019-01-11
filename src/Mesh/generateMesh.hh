//------------------------------------------------------------------------------
// Generate a mesh for the given set of NodeLists.
//------------------------------------------------------------------------------
#ifndef __Spheral_generateMesh__
#define __Spheral_generateMesh__

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension> class Mesh;

template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
void
generateMesh(const NodeListIterator nodeListBegin,
             const NodeListIterator nodeListEnd,
             const BoundaryIterator boundaryBegin,
             const BoundaryIterator boundaryEnd,
             const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const bool meshGhostNodes,
             const bool generateVoid,
             const bool generateParallelConnectivity,
             const bool removeBoundaryZones,
             const double voidThreshold,
             Mesh<Dimension>& mesh,
             NodeList<Dimension>& voidNodes);

}

#endif
