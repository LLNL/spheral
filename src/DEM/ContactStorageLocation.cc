//---------------------------------Spheral++----------------------------------//
// ContactIndex -- Simple structure used to track where pairwise variables 
//                 are stored in the pairFieldLists. 
//----------------------------------------------------------------------------//

#include <vector>

#include "DEM/ContactStorageLocation.hh"

namespace Spheral {

ContactIndex::ContactIndex():
  storeNodeList(0),
  storeNode(0),
  storeContact(0),
  pairNodeList(0),
  pairNode(0){
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
  pairNode(pairNodeIndex){
}


}