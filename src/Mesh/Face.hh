//---------------------------------Spheral++----------------------------------//
// Mesh::Face
//
// Created by JMO, Thu Oct 14 11:10:37 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Mesh_Face__
#define __Spheral_Mesh_Face__

namespace Spheral {

template<typename Dimension>
class Mesh<Dimension>::Face {
  //--------------------------- Public Interface ---------------------------//
public:

  //---------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //---------------------------------------------------------------------------
  Face(const Mesh<Dimension>& mesh,
       const unsigned ID,
       const int zone1ID,
       const int zone2ID,
       const std::vector<unsigned>& edgeIDs);

  // The ID of this Face.
  unsigned ID() const;

  // Sizes.
  unsigned numNodes() const;
  unsigned numEdges() const;

  // The ID's of the node and edges of the face.
  const std::vector<unsigned>& nodeIDs() const;
  const std::vector<unsigned>& edgeIDs() const;

  // ID's of the zones that share this face.
  int zone1ID() const;
  int zone2ID() const;

  // Position of the face.
  Vector position() const;

  // Area of the face (faceted).
  double area() const;

  // Unit normal of the face (such that nodes are counter-clockwise around the face).
  Vector unitNormal() const;

  // Return the other zone sharing this face.
  int oppositeZoneID(const int zoneID) const;

  // Is the given point above, below, or coplanar with the facet?
  int compare(const Vector& point,
              const double tol = 1.0e-8) const;

  //--------------------------- Private Interface ---------------------------//
private:
  const Mesh<Dimension>* mMeshPtr;
  unsigned mID;
  int mZone1ID, mZone2ID;
  std::vector<unsigned> mNodeIDs, mEdgeIDs;
  Face();

  friend class Mesh<Dimension>;
};

}

#include "FaceInline.hh"

#endif
