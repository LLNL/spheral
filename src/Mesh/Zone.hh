//---------------------------------Spheral++----------------------------------//
// Mesh::Zone
//
// Created by JMO, Thu Oct 14 11:10:37 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Mesh_Zone__
#define __Spheral_Mesh_Zone__

namespace Spheral {

template<typename Dimension>
class Mesh<Dimension>::Zone {
  //--------------------------- Public Interface ---------------------------//
public:

  //---------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //---------------------------------------------------------------------------
  Zone(const Mesh<Dimension>& mesh,
       const unsigned ID,
       const std::vector<int>& faceIDs);

  // ID of the zone.
  unsigned ID() const;

  // Sizes.
  unsigned numNodes() const;
  unsigned numEdges() const;
  unsigned numFaces() const;

  // The ID's of the subelements that make up the zone.
  const std::vector<unsigned>& nodeIDs() const;
  const std::vector<unsigned>& edgeIDs() const;
  const std::vector<int>& faceIDs() const;

  // Position.
  Vector position() const;

  // Volume.
  double volume() const;

  // Return the convex hull of the zone.
  ConvexHull convexHull() const;

  //--------------------------- Private Interface ---------------------------//
private:
  const Mesh<Dimension>* mMeshPtr;
  unsigned mID;
  std::vector<unsigned> mNodeIDs, mEdgeIDs;
  std::vector<int> mFaceIDs;
  Zone();

  friend class Mesh<Dimension>;
};

}

#include "ZoneInline.hh"

#endif
