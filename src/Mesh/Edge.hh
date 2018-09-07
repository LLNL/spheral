//---------------------------------Spheral++----------------------------------//
// Mesh::Edge
//
// Created by JMO, Thu Oct 14 11:10:37 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Mesh_Edge__
#define __Spheral_Mesh_Edge__

namespace Spheral {

template<typename Dimension>
class Mesh<Dimension>::Edge {
  //--------------------------- Public Interface ---------------------------//
public:

  //---------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //---------------------------------------------------------------------------
  Edge(const Mesh<Dimension>& mesh,
       const unsigned ID,
       const unsigned node1ID,
       const unsigned node2ID);
  unsigned ID() const;
  unsigned node1ID() const;
  unsigned node2ID() const;
  const Node& node1() const;
  const Node& node2() const;
  Vector position() const;
  double length() const;

  //--------------------------- Private Interface ---------------------------//
private:
  const Mesh<Dimension>* mMeshPtr;
  unsigned mID, mNode1ID, mNode2ID;
  Edge();

  friend class Mesh<Dimension>;
};

}

#include "EdgeInline.hh"

#endif
