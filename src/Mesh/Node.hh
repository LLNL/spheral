//---------------------------------Spheral++----------------------------------//
// Mesh::Node
//
// Created by JMO, Thu Oct 14 11:10:37 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Mesh_Node__
#define __Spheral_Mesh_Node__

namespace Spheral {

template<typename Dimension>
class Mesh<Dimension>::Node {
  //--------------------------- Public Interface ---------------------------//
public:

  //---------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //---------------------------------------------------------------------------
  Node(const Mesh<Dimension>& mesh,
       const unsigned ID,
       const std::vector<unsigned>& zoneIDs);
  unsigned ID() const;
  Vector position() const;

  // ID's of the zones that touch this node.
  const std::vector<unsigned>& zoneIDs() const;

  //--------------------------- Private Interface ---------------------------//
private:
  const Mesh<Dimension>* mMeshPtr;
  unsigned mID;
  std::vector<unsigned> mZoneIDs;
  Node();

  friend class Mesh<Dimension>;
};

}

#include "NodeInline.hh"

#endif
