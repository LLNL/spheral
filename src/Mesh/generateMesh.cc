//------------------------------------------------------------------------------
// Generate a mesh for the given set of NodeLists.
//------------------------------------------------------------------------------
#include <algorithm>

#include "generateMesh.hh"
#include "computeGenerators.hh"
#include "Mesh.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/generateVoidNodes.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/timingUtilities.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using BoundarySpace::Boundary;
using NodeSpace::NodeList;
using FieldSpace::Field;

template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
void
generateMesh(const NodeListIterator nodeListBegin,
             const NodeListIterator nodeListEnd,
             const BoundaryIterator boundaryBegin,
             const BoundaryIterator boundaryEnd,
             const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const bool generateVoid,
             const bool generateParallelConnectivity,
             const bool removeBoundaryZones,
             const double voidThreshold,
             Mesh<Dimension>& mesh,
             NodeList<Dimension>& voidNodes) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  // The total number of NodeLists we're working on.
  const size_t numNodeLists = distance(nodeListBegin, nodeListEnd);
  const unsigned voidOffset = distance(nodeListBegin, find(nodeListBegin, nodeListEnd, &voidNodes));

  // Pre-conditions.
  VERIFY2(voidOffset == numNodeLists - 1,
          "You must ensure the void nodes are part of the node set (at the end).");
  VERIFY2(voidNodes.numInternalNodes() == 0,
          "Void nodes not empty on startup!");
  VERIFY2(not (generateVoid and removeBoundaryZones),
          "You cannot simultaneously request generateVoid and removeBoundaryZones.");

  // Are we generating void?
  Timing::Time t0 = Timing::currentTime();
  vector<Vector> generators;
  vector<SymTensor> Hs;
  vector<unsigned> offsets;
  if (generateVoid or removeBoundaryZones) {
    // if (Process::getRank() == 0)
 cerr << "Computing void nodes" << endl;
    computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                     boundaryBegin, boundaryEnd,
                                                                     xmin, xmax, generators, Hs, offsets);
    unsigned numInternal = 0;
    double nPerh;
    for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd - 1; ++itr) {
      numInternal += (**itr).numInternalNodes();
      nPerh = (**itr).nodesPerSmoothingScale();
    }
    mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
    NodeSpace::generateVoidNodes(generators, Hs, mesh, xmin, xmax, numInternal, nPerh, voidThreshold, voidNodes);
  }

  // Extract the set of generators this domain needs (including any parallel neighbors).
  // This method gives us both the positions and Hs for the generators.
  // if (Process::getRank() == 0) 
cerr << "Computing generators" << endl;
  computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                   boundaryBegin, boundaryEnd,
                                                                   xmin, xmax, generators, Hs, offsets);
  MPI_Barrier(MPI_COMM_WORLD);
//   if (Process::getRank() == 0) 
cerr << "generateMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct generators." << endl;

  // Construct the mesh.
  t0 = Timing::currentTime();
  mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
  CHECK(mesh.numZones() == generators.size());
//   if (Process::getRank() == 0)
 cerr << "generateMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct mesh." << endl;

  // Remove any zones for generators that are not local to this domain.
  // if (Process::getRank() == 0) cerr << "Removing zones" << endl;
  t0 = Timing::currentTime();
  vector<unsigned> mask(mesh.numZones(), 0);
  unsigned ioff = 0;
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr, ++ioff) {
    if (ioff != voidOffset) {
      const unsigned zoneOffset = offsets[ioff];
      fill(mask.begin() + zoneOffset, mask.begin() + zoneOffset + (*itr)->numInternalNodes(), 1);
    }
  }
  if (not removeBoundaryZones) {
    fill(mask.begin() + offsets[voidOffset], mask.begin() + offsets[voidOffset] + voidNodes.numInternalNodes(), 1);
  }
  mesh.removeZonesByMask(mask);

  // If we removed boundary nodes, then the void was created for this purpose
  // alone.
  if (removeBoundaryZones) {
    voidNodes.numInternalNodes(0);
    offsets.back() = mesh.numZones();
  }
//   if (Process::getRank() == 0) 
cerr << "generateMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to remove boundary elements." << endl;

  // If requested we also compute the parallel connectivity.
  // if (Process::getRank() == 0) cerr << "Computing parallel connectivity" << endl;
  if (generateParallelConnectivity) {
    t0 = Timing::currentTime();
    mesh.generateDomainInfo();
//     if (Process::getRank() == 0) 
cerr << "generateMesh:: required " 
                                      << Timing::difference(t0, Timing::currentTime())
                                      << " seconds to generate parallel connectivity." << endl;
  }

  // Fill in the offset information.
  mesh.storeNodeListOffsets(nodeListBegin, nodeListEnd, offsets);

  // That's it.
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template void generateMesh<Dim<1>, 
                           vector<NodeList<Dim<1> >*>::iterator,
                           vector<Boundary<Dim<1> >*>::iterator>
              (const vector<NodeList<Dim<1> >*>::iterator nodeListBegin,
               const vector<NodeList<Dim<1> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
               const Dim<1>::Vector& xmin,
               const Dim<1>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<1> >& mesh,
               NodeList<Dim<1> >& voidNodes);

template void generateMesh<Dim<2>, 
                           vector<NodeList<Dim<2> >*>::iterator,
                           vector<Boundary<Dim<2> >*>::iterator>
              (const vector<NodeList<Dim<2> >*>::iterator nodeListBegin,
               const vector<NodeList<Dim<2> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
               const Dim<2>::Vector& xmin,
               const Dim<2>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<2> >& mesh,
               NodeList<Dim<2> >& voidNodes);

template void generateMesh<Dim<3>, 
                           vector<NodeList<Dim<3> >*>::iterator,
                           vector<Boundary<Dim<3> >*>::iterator>
              (const vector<NodeList<Dim<3> >*>::iterator nodeListBegin,
               const vector<NodeList<Dim<3> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
               const Dim<3>::Vector& xmin,
               const Dim<3>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<3> >& mesh,
               NodeList<Dim<3> >& voidNodes);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template void generateMesh<Dim<1>, 
                           vector<const NodeList<Dim<1> >*>::iterator,
                           vector<Boundary<Dim<1> >*>::const_iterator>
              (const vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
               const vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
               const Dim<1>::Vector& xmin,
               const Dim<1>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<1> >& mesh,
               NodeList<Dim<1> >& voidNodes);

template void generateMesh<Dim<2>, 
                           vector<const NodeList<Dim<2> >*>::iterator,
                           vector<Boundary<Dim<2> >*>::const_iterator>
              (const vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
               const vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
               const Dim<2>::Vector& xmin,
               const Dim<2>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<2> >& mesh,
               NodeList<Dim<2> >& voidNodes);

template void generateMesh<Dim<3>, 
                           vector<const NodeList<Dim<3> >*>::iterator,
                           vector<Boundary<Dim<3> >*>::const_iterator>
              (const vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
               const vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
               const Dim<3>::Vector& xmin,
               const Dim<3>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<3> >& mesh,
               NodeList<Dim<3> >& voidNodes);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template void generateMesh<Dim<1>, 
                           vector<const NodeList<Dim<1> >*>::iterator,
                           vector<Boundary<Dim<1> >*>::iterator>
              (const vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
               const Dim<1>::Vector& xmin,
               const Dim<1>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<1> >& mesh,
               NodeList<Dim<1> >& voidNodes);

template void generateMesh<Dim<2>, 
                           vector<const NodeList<Dim<2> >*>::iterator,
                           vector<Boundary<Dim<2> >*>::iterator>
              (const vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
               const Dim<2>::Vector& xmin,
               const Dim<2>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<2> >& mesh,
               NodeList<Dim<2> >& voidNodes);

template void generateMesh<Dim<3>, 
                           vector<const NodeList<Dim<3> >*>::iterator,
                           vector<Boundary<Dim<3> >*>::iterator>
              (const vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
               const vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
               const vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
               const vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
               const Dim<3>::Vector& xmin,
               const Dim<3>::Vector& xmax,
               const bool generateVoid,
               const bool generateParallelConnectivity,
               const bool removeBoundaryZones,
               const double voidThreshold,
               Mesh<Dim<3> >& mesh,
               NodeList<Dim<3> >& voidNodes);

}
}
