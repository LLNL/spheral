//------------------------------------------------------------------------------
// generateVoidNodes
//
// This algorithm tries to analyze how continuous a node distribution is, and 
// if it determines there is an edge to the distribution creates new void nodes
// outside that surface.
// We assume here that the caller has already created all the boundary ghost 
// nodes.
//------------------------------------------------------------------------------
#include <algorithm>
#include "boost/foreach.hpp"

#include "generateVoidNodes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/BSplineKernel.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Mesh/Mesh.hh"

namespace Spheral {
namespace NodeSpace {

using namespace std;
using std::abs;
using std::min;
using std::max;
using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using FieldSpace::FieldListBase;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using KernelSpace::BSplineKernel;
using MeshSpace::Mesh;

//------------------------------------------------------------------------------
// Overloaded method to find the closest point on a face for each dimension.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
closestPointOnFace(const Mesh<Dim<1> >& mesh,
                   const Mesh<Dim<1> >::Face& face,
                   const Dim<1>::Vector ri) {
  return face.position();
}

inline
Dim<2>::Vector
closestPointOnFace(const Mesh<Dim<2> >& mesh,
                   const Mesh<Dim<2> >::Face& face,
                   const Dim<2>::Vector ri) {
  const vector<unsigned>& nodeIDs = face.nodeIDs();
  REQUIRE(nodeIDs.size() == 2);
  return closestPointOnSegment(ri, 
                               mesh.node(nodeIDs[0]).position(),
                               mesh.node(nodeIDs[1]).position());
}

inline
Dim<3>::Vector
closestPointOnFace(const Mesh<Dim<3> >& mesh,
                   const Mesh<Dim<3> >::Face& face,
                   const Dim<3>::Vector ri) {
  return closestPointOnPlane(ri, 
                             mesh.node(face.nodeIDs()[0]).position(),
                             face.unitNormal());
}

//------------------------------------------------------------------------------
// The actual method.
//------------------------------------------------------------------------------
template<typename Dimension>
void generateVoidNodes(const vector<typename Dimension::Vector>& generators,
                       const vector<typename Dimension::SymTensor>& Hs,
                       const Mesh<Dimension>& mesh,
                       const typename Dimension::Vector& xmin,
                       const typename Dimension::Vector& xmax,
                       const unsigned numInternal,
                       const double nPerh,
                       const double threshold,
                       NodeList<Dimension>& voidNodes) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ConvexHull ConvexHull;
  typedef typename Mesh<Dimension>::Face Face;

  const double tiny = 1.0e-5;

  // Pre-conditions.
  VERIFY(generators.size() >= numInternal);
  VERIFY(generators.size() == Hs.size());
  VERIFY(mesh.numZones() >= numInternal);
  VERIFY(mesh.numZones() <= generators.size());
  VERIFY(voidNodes.numNodes() == 0);
  VERIFY(nPerh > 0.0);
  VERIFY(threshold >= 0.0);

  unsigned i, j, k, ii, ivoid = 0;
  Vector ri, rj, rij, etai, etaj, fmhat;
  SymTensor Hi, Hj;
  Field<Dimension, Vector>& vpos = voidNodes.positions();
  Field<Dimension, SymTensor>& vH = voidNodes.Hfield();

  // For each face, we build the generator positions and H's on each side of the face.
  vector<pair<Vector, Vector> > faceGenerators(mesh.numFaces(),
                                               make_pair(Vector::one * 1e100,
                                                         Vector::one * 1e100));
  vector<pair<SymTensor, SymTensor> > faceHs(mesh.numFaces(),
                                             make_pair(SymTensor::one,
                                                       SymTensor::one));
  for (i = 0; i != mesh.numFaces(); ++i) {
    const Face& face = mesh.face(i);
    CHECK2(Mesh<Dimension>::positiveID(face.zone1ID()) == Mesh<Dimension>::UNSETID or Mesh<Dimension>::positiveID(face.zone1ID()) < generators.size(), i << " " << face.zone1ID());
    CHECK2(Mesh<Dimension>::positiveID(face.zone2ID()) == Mesh<Dimension>::UNSETID or Mesh<Dimension>::positiveID(face.zone2ID()) < generators.size(), i << " " << face.zone2ID());
    if (Mesh<Dimension>::positiveID(face.zone1ID()) != Mesh<Dimension>::UNSETID) {
      j = Mesh<Dimension>::positiveID(face.zone1ID());
      faceGenerators[i].first = generators[j];
      faceHs[i].first = Hs[j];
    }
    if (Mesh<Dimension>::positiveID(face.zone2ID()) != Mesh<Dimension>::UNSETID) {
      j = Mesh<Dimension>::positiveID(face.zone2ID());
      faceGenerators[i].second = generators[j];
      faceHs[i].second = Hs[j];
    }
  }

#ifdef USE_MPI
  // In parallel we need the other domains generator positions and H's for shared faces.
  const vector<unsigned>& neighborDomains = mesh.neighborDomains();
  const vector<vector<unsigned> >& sharedFaces = mesh.sharedFaces();
  const unsigned numNeighbors = neighborDomains.size();
  CHECK(sharedFaces.size() == numNeighbors);

  // Post all our sends.
  vector<unsigned> bufSizes;
  list<vector<char> > sendBufs;
  vector<MPI_Request> sendRequests;
  sendRequests.reserve(numNeighbors);
  for (unsigned irank = 0; irank != numNeighbors; ++irank) {
    if (sharedFaces[irank].size() > 0) {
      sendBufs.push_back(vector<char>());
      for (typename vector<unsigned>::const_iterator itr = sharedFaces[irank].begin();
           itr != sharedFaces[irank].end();
           ++itr) {
        const Face& face = mesh.face(*itr);
        CHECK(Mesh<Dimension>::positiveID(face.zone1ID()) == Mesh<Dimension>::UNSETID or
              Mesh<Dimension>::positiveID(face.zone2ID()) == Mesh<Dimension>::UNSETID);
        i = (Mesh<Dimension>::positiveID(face.zone1ID()) == Mesh<Dimension>::UNSETID ? 
             Mesh<Dimension>::positiveID(face.zone2ID()) :
             Mesh<Dimension>::positiveID(face.zone1ID()));
        CHECK(i < generators.size());
        packElement(generators[i], sendBufs.back());
        packElement(Hs[i], sendBufs.back());
      }
      bufSizes.push_back(sendBufs.back().size());
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&sendBufs.back().front(), bufSizes.back(), MPI_CHAR, neighborDomains[irank], 1, MPI_COMM_WORLD, &sendRequests.back());
    }
  }

  // Now go through and get all our receives from our neighbors.
  unsigned ioff = 0;
  for (unsigned irank = 0; irank != numNeighbors; ++irank) {
    if (sharedFaces[irank].size() > 0) {
      vector<char> buffer(bufSizes[ioff]);
      MPI_Status recvStatus;
      MPI_Recv(&buffer.front(), bufSizes[ioff], MPI_CHAR, neighborDomains[irank], 1, MPI_COMM_WORLD, &recvStatus);
      vector<char>::const_iterator bufItr = buffer.begin();
      for (typename vector<unsigned>::const_iterator itr = sharedFaces[irank].begin();
           itr != sharedFaces[irank].end();
           ++itr) {
        const Face& face = mesh.face(*itr);
        unpackElement(rj, bufItr, buffer.end());
        unpackElement(Hj, bufItr, buffer.end());
        if (Mesh<Dimension>::positiveID(face.zone1ID()) == Mesh<Dimension>::UNSETID) {
          faceGenerators[*itr].first = rj;
          faceHs[*itr].first = Hj;
        } else {
          faceGenerators[*itr].second = rj;
          faceHs[*itr].second = Hj;
        }
      }        
      CHECK(bufItr == buffer.end());
    }
  }
#endif

  // Now we can walk the generators.
  for (i = 0; i != numInternal; ++i) {
    ri = generators[i];
    Hi = Hs[i];

    // // Generate the convex hull of this zones and it's neighbors.
    // const vector<unsigned>& nodeIDs = mesh.zone(i).nodeIDs();
    // vector<Vector> neighborNodes;
    // for (k = 0; k != nodeIDs.size(); ++k) {
    //   const vector<unsigned>& zoneIDs = mesh.node(nodeIDs[k]).zoneIDs();
    //   for (j = 0; j != zoneIDs.size(); ++j) {
    //     CHECK(zoneIDs[j] < mesh.numZones());
    //     const vector<unsigned>& otherNodeIDs = mesh.zone(zoneIDs[j]).nodeIDs();
    //     for (ii = 0; ii != otherNodeIDs.size(); ++ii) {
    //       neighborNodes.push_back(mesh.node(otherNodeIDs[ii]).position());
    //     }
    //   }
    // }
    // const ConvexHull neighborHull(neighborNodes);

    // Look for any faces of this node's cell across which we 
    // should be generating a void point.
    const vector<int>& faceIDs = mesh.zone(i).faceIDs();
    BOOST_FOREACH(int faceID, faceIDs) {
      j = Mesh<Dimension>::positiveID(faceID);
      const typename Mesh<Dimension>::Face& face = mesh.face(j);
      rj = (Mesh<Dimension>::positiveID(face.zone1ID()) == i ? faceGenerators[j].second : faceGenerators[j].first);
      Hj = (Mesh<Dimension>::positiveID(face.zone1ID()) == i ? faceHs[j].second : faceHs[j].first);
      rij = ri - rj;
      etai = Hi*rij;
      etaj = Hj*rij;
      if (etai.magnitude() > threshold and etaj.magnitude() > threshold) {

        // Yep, generate a new void point.
        fmhat = (face.position() - ri).unitVector(); //  (closestPointOnFace(mesh, face, ri) - ri).unitVector();
        rj = ri + Hi.Inverse()/nPerh * fmhat;
        if (testPointInBox(rj, xmin, xmax)) {
          voidNodes.numInternalNodes(ivoid + 1);
          vpos[ivoid] = rj;
          vH[ivoid] = Hi;
          ++ivoid;
        }
      }
    }
  }

#ifdef USE_MPI
  // Wait until our sends are all satisfied.
  vector<MPI_Status> sendStatus(sendRequests.size());
  MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
#endif

  // That's it.
}

}
}
