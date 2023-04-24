//---------------------------------Spheral++----------------------------------//
// ContactIndex -- Simple structure used to track where pairwise variables 
//                 are stored in the pairFieldLists. A -1 indicates an invalid
//                 entry. The contact index has two configurations. The first
//                 is for particle-particle contacts, for which store* and 
//                 pair* properties are used. The second is for particle-
//                 solid boundary interactions and there store* and 
//                 solidBoundary are used. In both case, unused properties are
//                 flagged as invalid (-1).
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "DEM/ContactStorageLocation.hh"

namespace Spheral {

ContactIndex::ContactIndex():
  storeNodeList(-1),
  storeNode(-1),
  storeContact(-1),
  pairNodeList(-1),
  pairNode(-1),
  solidBoundary(-1){
}

ContactIndex::ContactIndex(const int storeNodeListIndex, 
                           const int storeNodeIndex, 
                           const int storeContactIndex,
                           const int pairNodeListIndex,
                           const int pairNodeIndex):
  storeNodeList(storeNodeListIndex),
  storeNode(storeNodeIndex),
  storeContact(storeContactIndex),
  pairNodeList(pairNodeListIndex),
  pairNode(pairNodeIndex),
  solidBoundary(-1){
}

ContactIndex::ContactIndex(const int storeNodeListIndex, 
                           const int storeNodeIndex, 
                           const int storeContactIndex,
                           const int solidBoundaryIndex):
  storeNodeList(storeNodeListIndex),
  storeNode(storeNodeIndex),
  storeContact(storeContactIndex),
  pairNodeList(-1),
  pairNode(-1),
  solidBoundary(solidBoundaryIndex){
}



}