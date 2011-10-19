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

  // Look for any boundaries whose ghost nodes we should include.
  for (BoundaryIterator bcItr = boundaryBegin;
       bcItr != boundaryEnd;
       ++bcItr) {
    for (NodeListIterator nodeListItr = nodeListBegin; nodeListItr != nodeListEnd; ++nodeListItr) {
      if ((*bcItr)->meshGhostNodes() and (*bcItr)->haveNodeList(**nodeListItr)) {
        const Field<Dimension, Vector>& pos = (**nodeListItr).positions();
        const Field<Dimension, SymTensor>& H = (**nodeListItr).Hfield();
        const vector<int>& ghostNodes = (*bcItr)->ghostNodes(**nodeListItr);
        for (k = 0; k != ghostNodes.size(); ++k) {
          i = ghostNodes[k];
//           if (!(i < pos.numElements())) {
//             cerr << "Blago!  " << ghostNodes.size() << endl << "  ---> ";
//             for (unsigned kk = 0; kk != ghostNodes.size(); ++kk) cerr << ghostNodes[kk] << " ";
//             cerr << endl;
//           }
          CHECK2(i < pos.numElements(), "Out of bounds:  " << i << " " << pos.numElements());
          if (testPointInBox(pos(i), xmin, xmax)) {
            ++offsets.back();
            ++nlocal;
            localPositions.push_back(pos(i));
            localHs.push_back(H(i));
          }
        }
      }
    }
  }
  CHECK(localPositions.size() == nlocal);
  CHECK(localHs.size() == nlocal);
  CHECK(offsets.size() == numNodeLists + 1);

  // Copy our local info to the result arrays.
  positions = localPositions;
  Hs = localHs;

#ifdef USE_MPI
  if (numDomains > 1) {

    // Compute the convex hull of each domain, and distribute them to everyone.
//     if (Process::getRank() == 0) cerr << "Computing and broadcasting local hulls of domains." << endl;
    const ConvexHull localHull(localPositions);
    vector<ConvexHull> domainHulls(numDomains, localHull);
    vector<unsigned> domainZoneOffset(1, 0);
    {
      vector<char> localBuffer;
      packElement(localHull, localBuffer);
      for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
        vector<char> buffer = localBuffer;
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
        vector<char>::const_iterator itr = buffer.begin();
        unpackElement(domainHulls[sendProc], itr, buffer.end());
        CHECK(itr == buffer.end());
        domainZoneOffset.push_back(domainZoneOffset.back() + domainHulls[sendProc].vertices().size());
      }
    }
    CHECK(domainHulls.size() == numDomains);
    CHECK(domainZoneOffset.size() == numDomains + 1);

    // Create a mesh of the hull points for all domains.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Computing mesh of union of local hulls." << endl;
    vector<Vector> hullGenerators;
    for (unsigned k = 0; k != domainHulls.size(); ++k) {
      const vector<Vector>& hullVertices = domainHulls[k].vertices();
      copy(hullVertices.begin(), hullVertices.end(), back_inserter(hullGenerators));
    }
    Mesh<Dimension> hullMesh(hullGenerators, xmin, xmax);

    // Build up the set of domains we need to communicate with according to two criteria:
    //  1.  Any domain hull that intersects our own.
    //  2.  Any domain hull that has elements adjacent to one of ours in the hullMesh.
    set<unsigned> neighborSet;

    // First any hulls that intersect ours.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Checking for hulls that intersect local." << endl;
    for (unsigned otherProc = 0; otherProc != numDomains; ++otherProc) {
      if (otherProc != rank and
          localHull.intersect(domainHulls[otherProc])) neighborSet.insert(otherProc);
    }

    // Now any hulls that have elements adjacent to ours in the hull mesh.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Checking for hulls adjacent in mesh." << endl;
    for (unsigned izone = domainZoneOffset[rank];
         izone != domainZoneOffset[rank + 1];
         ++izone) {
      const vector<unsigned>& nodeIDs = hullMesh.zone(izone).nodeIDs();
      for (typename vector<unsigned>::const_iterator nodeItr = nodeIDs.begin();
           nodeItr != nodeIDs.end();
           ++nodeItr) {
        const unsigned inode = *nodeItr;
        const vector<unsigned>& nodeZoneIDs = hullMesh.node(inode).zoneIDs();
        for (typename vector<unsigned>::const_iterator zoneItr = nodeZoneIDs.begin();
             zoneItr != nodeZoneIDs.end();
             ++zoneItr) {
          const unsigned izoneNeighbor = *zoneItr;
          if (izoneNeighbor != Mesh<Dimension>::UNSETID and
              (izoneNeighbor < domainZoneOffset[rank] or izoneNeighbor >= domainZoneOffset[rank + 1])) {
            const int otherProc = bisectSearch(domainZoneOffset, izoneNeighbor);
            CHECK(otherProc >= 0 and otherProc < domainZoneOffset.size() - 1);
            CHECK(izoneNeighbor >= domainZoneOffset[otherProc] and
                  izoneNeighbor <  domainZoneOffset[otherProc + 1]);
            neighborSet.insert(unsigned(otherProc));
          }
        }
      }
    }

    // Make sure everyone is consistent about who talks to whom.
    vector<unsigned> neighborDomains;
    copy(neighborSet.begin(), neighborSet.end(), back_inserter(neighborDomains));
    sort(neighborDomains.begin(), neighborDomains.end());
    BEGIN_CONTRACT_SCOPE;
    {
      for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
        unsigned numOthers = neighborDomains.size();
        vector<unsigned> otherNeighbors(neighborDomains);
        MPI_Bcast(&numOthers, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (numOthers > 0) {
          otherNeighbors.resize(numOthers);
          MPI_Bcast(&(otherNeighbors.front()), numOthers, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
          CHECK2(rank == sendProc or
                 count(neighborDomains.begin(), neighborDomains.end(), sendProc) == 
                 count(otherNeighbors.begin(), otherNeighbors.end(), rank),
                 "Bad neighbor connectivity:  "
                 << count(neighborDomains.begin(), neighborDomains.end(), sendProc) << " != "
                 << count(otherNeighbors.begin(), otherNeighbors.end(), rank));
        }
      }
    }
    END_CONTRACT_SCOPE;

    // Build the mesh of the local generators, and extract the convex hulls of each zone.
    const Mesh<Dimension> localMesh(localPositions, xmin, xmax);
    CHECK(localMesh.numZones() == nlocal);
    vector<ConvexHull> localZoneHulls;
    for (unsigned i = 0; i != nlocal; ++i) localZoneHulls.push_back(localMesh.zone(i).convexHull());

    // Now we have to determine which of our generators go to each neighbor.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Sending local generators to neighbors." << endl;
    list<vector<char> > localBuffers;
    list<unsigned> localBufSizes;
    vector<MPI_Request> sendRequests(2*neighborDomains.size());
    for (unsigned k = 0; k != neighborDomains.size(); ++k) {
      const unsigned otherProc = neighborDomains[k];

      // Look for any zones that intersect the other domain hull.  We'll send that generator
      // along with its immediate neighbors.
      set<unsigned> sendGenIDs;
      unsigned j;
      for (unsigned i = 0; i != nlocal; ++i) {
        if (domainHulls[otherProc].convexIntersect(localZoneHulls[i])) {
          sendGenIDs.insert(i);
          const vector<unsigned>& faceIDs = localMesh.zone(i).faceIDs();
          for (vector<unsigned>::const_iterator faceItr = faceIDs.begin();
               faceItr != faceIDs.end();
               ++faceItr) {
            j = localMesh.face(*faceItr).oppositeZoneID(i);
            if (j != Mesh<Dimension>::UNSETID) sendGenIDs.insert(j);
          }
        }
      }
//       if (Process::getRank() == 0) cerr << "    computeGenerators:  " << Process::getRank() << "->" << otherProc << " sending " << sendGenIDs.size() << " of " << localPositions.size() << endl;

      // Pack up the local generators we're sending.
      vector<Vector> sendgens;
      vector<SymTensor> sendHs;
      for (typename set<unsigned>::const_iterator itr = sendGenIDs.begin();
           itr != sendGenIDs.end();
           ++itr) {
        sendgens.push_back(localPositions[*itr]);
        sendHs.push_back(localHs[*itr]);
      }
      CHECK(sendgens.size() == sendGenIDs.size());
      CHECK(sendHs.size() == sendGenIDs.size());

//       // Blago!
//       for (unsigned i = 0; i != localPositions.size(); ++i) {
//         sendgens.push_back(localPositions[i]);
//         sendHs.push_back(localHs[i]);
//       }
//       // Blago!

      localBuffers.push_back(vector<char>());
      packElement(sendgens, localBuffers.back());
      packElement(sendHs, localBuffers.back());
      localBufSizes.push_back(localBuffers.back().size());

      // Send the generator info.
      MPI_Isend(&localBufSizes.back(), 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &(sendRequests[2*k]));
      MPI_Isend(&localBuffers.back().front(), localBufSizes.back(), MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &(sendRequests[2*k + 1]));
    }

    // Get the info from each of our neighbors and append it to the result.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Receiving neighbor generators." << endl;
    for (unsigned k = 0; k != neighborDomains.size(); ++k) {
      const unsigned recvProc = neighborDomains[k];
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, recvProc, 1, MPI_COMM_WORLD, &recvStatus1);
      CHECK(bufSize > 0);
      vector<char> buffer(bufSize, '\0');
      MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, recvProc, 2, MPI_COMM_WORLD, &recvStatus2);
      vector<char>::const_iterator itr = buffer.begin();
      unpackElement(positions, itr, buffer.end());
      unpackElement(Hs, itr, buffer.end());
      CHECK(itr == buffer.end());
    }

    // Make sure all our sends are completed.
//     if (Process::getRank() == 0) cerr << "computeGenerators:  Done." << endl;
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
  }
#endif

  // That's it.
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template void computeGenerators<Dim<1>, 
                                vector<NodeList<Dim<1> >*>::iterator,
                                vector<Boundary<Dim<1> >*>::iterator>
         (vector<NodeList<Dim<1> >*>::iterator nodeListBegin,
          vector<NodeList<Dim<1> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
          const Dim<1>::Vector& xmin,
          const Dim<1>::Vector& xmax,
          vector<Dim<1>::Vector>& positions,
          vector<Dim<1>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<2>, 
                                vector<NodeList<Dim<2> >*>::iterator,
                                vector<Boundary<Dim<2> >*>::iterator>
         (vector<NodeList<Dim<2> >*>::iterator nodeListBegin,
          vector<NodeList<Dim<2> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
          const Dim<2>::Vector& xmin,
          const Dim<2>::Vector& xmax,
          vector<Dim<2>::Vector>& positions,
          vector<Dim<2>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<3>, 
                                vector<NodeList<Dim<3> >*>::iterator,
                                vector<Boundary<Dim<3> >*>::iterator>
         (vector<NodeList<Dim<3> >*>::iterator nodeListBegin,
          vector<NodeList<Dim<3> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
          const Dim<3>::Vector& xmin,
          const Dim<3>::Vector& xmax,
          vector<Dim<3>::Vector>& positions,
          vector<Dim<3>::SymTensor>& Hs,
          vector<unsigned>& offsets);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template void computeGenerators<Dim<1>, 
                                vector<const NodeList<Dim<1> >*>::iterator,
                                vector<Boundary<Dim<1> >*>::iterator>
         (vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
          const Dim<1>::Vector& xmin,
          const Dim<1>::Vector& xmax,
          vector<Dim<1>::Vector>& positions,
          vector<Dim<1>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<2>, 
                                vector<const NodeList<Dim<2> >*>::iterator,
                                vector<Boundary<Dim<2> >*>::iterator>
         (vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
          const Dim<2>::Vector& xmin,
          const Dim<2>::Vector& xmax,
          vector<Dim<2>::Vector>& positions,
          vector<Dim<2>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<3>, 
                                vector<const NodeList<Dim<3> >*>::iterator,
                                vector<Boundary<Dim<3> >*>::iterator>
         (vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
          vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
          const Dim<3>::Vector& xmin,
          const Dim<3>::Vector& xmax,
          vector<Dim<3>::Vector>& positions,
          vector<Dim<3>::SymTensor>& Hs,
          vector<unsigned>& offsets);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template void computeGenerators<Dim<1>, 
                                vector<const NodeList<Dim<1> >*>::iterator,
                                vector<Boundary<Dim<1> >*>::const_iterator>
         (vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
          const Dim<1>::Vector& xmin,
          const Dim<1>::Vector& xmax,
          vector<Dim<1>::Vector>& positions,
          vector<Dim<1>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<2>, 
                                vector<const NodeList<Dim<2> >*>::iterator,
                                vector<Boundary<Dim<2> >*>::const_iterator>
         (vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
          const Dim<2>::Vector& xmin,
          const Dim<2>::Vector& xmax,
          vector<Dim<2>::Vector>& positions,
          vector<Dim<2>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<3>, 
                                vector<const NodeList<Dim<3> >*>::iterator,
                                vector<Boundary<Dim<3> >*>::const_iterator>
         (vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
          vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
          vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
          const Dim<3>::Vector& xmin,
          const Dim<3>::Vector& xmax,
          vector<Dim<3>::Vector>& positions,
          vector<Dim<3>::SymTensor>& Hs,
          vector<unsigned>& offsets);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template void computeGenerators<Dim<1>, 
                                vector<NodeList<Dim<1> >*>::const_iterator,
                                vector<Boundary<Dim<1> >*>::const_iterator>
         (vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
          vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
          vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
          const Dim<1>::Vector& xmin,
          const Dim<1>::Vector& xmax,
          vector<Dim<1>::Vector>& positions,
          vector<Dim<1>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<2>, 
                                vector<NodeList<Dim<2> >*>::const_iterator,
                                vector<Boundary<Dim<2> >*>::const_iterator>
         (vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
          vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
          vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
          const Dim<2>::Vector& xmin,
          const Dim<2>::Vector& xmax,
          vector<Dim<2>::Vector>& positions,
          vector<Dim<2>::SymTensor>& Hs,
          vector<unsigned>& offsets);

template void computeGenerators<Dim<3>, 
                                vector<NodeList<Dim<3> >*>::const_iterator,
                                vector<Boundary<Dim<3> >*>::const_iterator>
         (vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
          vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
          vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
          vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
          const Dim<3>::Vector& xmin,
          const Dim<3>::Vector& xmax,
          vector<Dim<3>::Vector>& positions,
          vector<Dim<3>::SymTensor>& Hs,
          vector<unsigned>& offsets);

}
}

