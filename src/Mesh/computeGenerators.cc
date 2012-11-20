//---------------------------------Spheral++----------------------------------//
// computeGenerators.
//
// Helper method used when we're generating meshes.
// This method takes a set of NodeLists, figures out what domains have to talk
// to one another, and returns the flattened set of positions and Hs for the 
// the generators this domain needs (including those from neighbor domains).
//
// Created by JMO, Mon Dec  6 10:34:43 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <set>

#include "boost/foreach.hpp"

#include "computeGenerators.hh"
#include "Mesh.hh"
#include "MeshConstructionUtilities.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/packElement.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/removeElements.hh"

#ifdef USE_MPI
extern "C" {
#include "mpi.h"
}
#endif

namespace Spheral {
namespace MeshSpace {

using namespace std;
using NodeSpace::NodeList;
using BoundarySpace::Boundary;
using FieldSpace::Field;
using NeighborSpace::Neighbor;

// //------------------------------------------------------------------------------
// // A helper to build up the hull generators.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// conditionalInsertPoint(set<typename Mesh<Dimension>::Key>& hullGenerators,
//                        const typename Dimension::Vector& x,
//                        const typename Dimension::Vector& xmin,
//                        const typename Dimension::Vector& xmax,
//                        const typename Dimension::Vector& boxInv) {
//   if (testPointInBox(x, xmin, xmax)) hullGenerators.insert(hashPosition(x, xmin, xmax, boxInv));
// }

// void
// addBoundingPoints(set<Mesh<Dim<1> >::Key>& hullGenerators,
//                   const Dim<1>::Vector& ri,
//                   const Dim<1>::Vector& extenti,
//                   const Dim<1>::Vector& xmin,
//                   const Dim<1>::Vector& xmax,
//                   const Dim<1>::Vector& boxInv) {
//   conditionalInsertPoint<Dim<1> >(hullGenerators, ri, xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<1> >(hullGenerators, ri - extenti, xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<1> >(hullGenerators, ri + extenti, xmin, xmax, boxInv);
// }

// void
// addBoundingPoints(set<Mesh<Dim<2> >::Key>& hullGenerators,
//                   const Dim<2>::Vector& ri,
//                   const Dim<2>::Vector& extenti,
//                   const Dim<2>::Vector& xmin,
//                   const Dim<2>::Vector& xmax,
//                   const Dim<2>::Vector& boxInv) {
//   typedef Dim<2>::Vector Vector;
//   conditionalInsertPoint<Dim<2> >(hullGenerators, ri, xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<2> >(hullGenerators, ri - extenti, xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<2> >(hullGenerators, ri + extenti, xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<2> >(hullGenerators, ri + Vector( extenti.x(), -extenti.y()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<2> >(hullGenerators, ri + Vector(-extenti.x(),  extenti.y()), xmin, xmax, boxInv);
// }

// void
// addBoundingPoints(set<Mesh<Dim<3> >::Key>& hullGenerators,
//                   const Dim<3>::Vector& ri,
//                   const Dim<3>::Vector& extenti,
//                   const Dim<3>::Vector& xmin,
//                   const Dim<3>::Vector& xmax,
//                   const Dim<3>::Vector& boxInv) {
//   typedef Dim<3>::Vector Vector;
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri,                                                    xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector( extenti.x(),  extenti.y(), -extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector(-extenti.x(),  extenti.y(), -extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector( extenti.x(), -extenti.y(), -extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector(-extenti.x(), -extenti.y(), -extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector( extenti.x(),  extenti.y(),  extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector(-extenti.x(),  extenti.y(),  extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector( extenti.x(), -extenti.y(),  extenti.z()), xmin, xmax, boxInv);
//   conditionalInsertPoint<Dim<3> >(hullGenerators, ri + Vector(-extenti.x(), -extenti.y(),  extenti.z()), xmin, xmax, boxInv);
// }

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
void
computeGenerators(NodeListIterator nodeListBegin,
                  NodeListIterator nodeListEnd,
                  BoundaryIterator boundaryBegin,
                  BoundaryIterator boundaryEnd,
                  const typename Dimension::Vector& xmin,
                  const typename Dimension::Vector& xmax,
                  const bool generateParallelRind,
                  vector<typename Dimension::Vector>& positions,
                  vector<typename Dimension::SymTensor>& Hs,
                  vector<unsigned>& offsets) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ConvexHull ConvexHull;
  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Key Key;

  // Parallel geometry.
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  const unsigned numNodeLists = distance(nodeListBegin, nodeListEnd);

  // Flatten the local positions and Hs.
  vector<Vector> localPositions;
  vector<SymTensor> localHs;
  offsets = vector<unsigned>(1, 0);
  unsigned i, k, nlocal = 0;
  for (NodeListIterator nodeListItr = nodeListBegin; nodeListItr != nodeListEnd; ++nodeListItr) {
    const Field<Dimension, Vector>& pos = (**nodeListItr).positions();
    const Field<Dimension, SymTensor>& H = (**nodeListItr).Hfield();
    copy(pos.internalBegin(), pos.internalEnd(), back_inserter(localPositions));
    copy(H.internalBegin(), H.internalEnd(), back_inserter(localHs));
    offsets.push_back(offsets.back() + (**nodeListItr).numInternalNodes());
    nlocal += (**nodeListItr).numInternalNodes();
  }

//   // Look for any boundaries whose ghost nodes we should include.
//   for (BoundaryIterator bcItr = boundaryBegin;
//        bcItr != boundaryEnd;
//        ++bcItr) {
//     for (NodeListIterator nodeListItr = nodeListBegin; nodeListItr != nodeListEnd; ++nodeListItr) {
//       if ((*bcItr)->meshGhostNodes() and (*bcItr)->haveNodeList(**nodeListItr)) {
//         const Field<Dimension, Vector>& pos = (**nodeListItr).positions();
//         const Field<Dimension, SymTensor>& H = (**nodeListItr).Hfield();
//         const vector<int>& ghostNodes = (*bcItr)->ghostNodes(**nodeListItr);
//         for (k = 0; k != ghostNodes.size(); ++k) {
//           i = ghostNodes[k];
// //           if (!(i < pos.numElements())) {
// //             cerr << "Blago!  " << ghostNodes.size() << endl << "  ---> ";
// //             for (unsigned kk = 0; kk != ghostNodes.size(); ++kk) cerr << ghostNodes[kk] << " ";
// //             cerr << endl;
// //           }
//           CHECK2(i < pos.numElements(), "Out of bounds:  " << i << " " << pos.numElements());
//           if (testPointInBox(pos(i), xmin, xmax)) {
//             ++offsets.back();
//             ++nlocal;
//             localPositions.push_back(pos(i));
//             localHs.push_back(H(i));
//           }
//         }
//       }
//     }
//   }
  CHECK(localPositions.size() == nlocal);
  CHECK(localHs.size() == nlocal);
  CHECK(offsets.size() == numNodeLists + 1);

  // Copy our local info to the result arrays.
  positions = localPositions;
  Hs = localHs;

#ifdef USE_MPI
  // If requested we can generate the parallel rind of generators.
  if (generateParallelRind and numDomains > 1) {

    // Build up the mesh of our local generators.
    Mesh<Dimension> localMesh(localPositions, xmin, xmax);
    const vector<unsigned>& neighborDomains = localMesh.neighborDomains();
    const vector<vector<unsigned> >& sharedNodes = localMesh.sharedNodes();
    const unsigned numNeighborDomains = neighborDomains.size();

    // Tell every domain we share a node with about our cells that have that node.
    list<vector<char> > sendBufs;
    vector<unsigned> sendSizes(numNeighborDomains);
    vector<MPI_Request> sendRequests(2*numNeighborDomains);
    for (unsigned kdomain = 0; kdomain != numNeighborDomains; ++kdomain) {
      const unsigned otherProc = neighborDomains[kdomain];
      CHECK(sharedNodes[kdomain].size() > 0);

      // Pack up our generators for the other domain.
      sendBufs.push_back(vector<char>());
      vector<char>& buf = sendBufs.back();
      for (vector<unsigned>::const_iterator nodeItr = sharedNodes[kdomain].begin();
           nodeItr != sharedNodes[kdomain].end();
           ++nodeItr) {
        const vector<unsigned>& cells = localMesh.node(*nodeItr).zoneIDs();
        for (vector<unsigned>::const_iterator cellItr = cells.begin();
             cellItr != cells.end();
             ++cellItr) {
          if (Mesh<Dimension>::positiveID(*cellItr) != Mesh<Dimension>::UNSETID) {
            CHECK2(*cellItr < localPositions.size(), *cellItr << " " << localPositions.size());
            packElement(localPositions[*cellItr], buf);
            packElement(localHs[*cellItr], buf);
          }
        }
      }
      sendSizes[kdomain] = buf.size();
      MPI_Isend(&sendSizes[kdomain], 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &sendRequests[2*kdomain]);
      MPI_Isend(&buf.front(), buf.size(), MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &sendRequests[2*kdomain+1]);
    }
    CHECK(sendBufs.size() == numNeighborDomains);

    // Gather up the neighbor generators for each node we share with them.
    // We rely upon the fact that these generators will be unique for this step!
    for (unsigned kdomain = 0; kdomain != numNeighborDomains; ++kdomain) {
      const unsigned otherProc = neighborDomains[kdomain];
      CHECK(sharedNodes[kdomain].size() > 0);
      MPI_Status status1, status2;
      unsigned bufSize;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &status1);
      CHECK(bufSize > 0);
      vector<char> buffer(bufSize);
      MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &status2);
      vector<char>::const_iterator bufItr = buffer.begin();
      Vector xi;
      SymTensor Hi;
      while (bufItr != buffer.end()) {
        unpackElement(xi, bufItr, buffer.end());
        unpackElement(Hi, bufItr, buffer.end());
        positions.push_back(xi);
        Hs.push_back(Hi);
      }
    }

    // Wait until all our sends have been satisfied.
    vector<MPI_Status> status(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &status.front());
  }
#endif

  // That's it.
}

}
}

