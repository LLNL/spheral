#include <vector>

#include "DEM/ContactStorageLocation.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// constructor for our contact index
//------------------------------------------------------------------------------
ContactIndex::ContactIndex(const int nodeListIndex, 
                           const int nodeIndex, 
                           const int contactIndex):
  nodeListIndex(nodeListIndex),
  nodeIndex(nodeIndex),
  contactIndex(contactIndex){

}

//------------------------------------------------------------------------------
// constructor for our contact index
//------------------------------------------------------------------------------
ContactStorageLocation::ContactStorageLocation(const int nodeListi,
                                               const int i,
                                               const int uniqueIndexi,
                                               const int numInternalNodesi,
                                               const std::vector<int>& contactsi,
                                               const int nodeListj,
                                               const int j,
                                               const int uniqueIndexj,
                                               const int numInternalNodesj,
                                               const std::vector<int>& contactsj)
  storageLocation(),
  pairLocation(){

    const auto selectNodei = (i < numInternalNodesi) and (uniqueIndexi < uniqueIndexj or j > numInternalNodesj);
    int storageNodeListIndex, storageNodeIndex, uniqueIndex;
    if (selectNodei) {
      storageNodeListIndex = nodeListi;
      storageNodeIndex = i;
      uniqueIndex = uniqueIndexj;
    } else{
      storageNodeListIndex = nodeListj;
      storageNodeIndex = j;
      uniqueIndex = uniqueIndexi;
    }

    // return the nodelist,node ids plus the unique index we'll search for
    std::vector<int> storageIndices = {storageNodeListIndex,storageNodeIndex,uniqueIndex};
    auto& neighborContacts = neighborIndices(storageNodeListIndex,storageNodeIndex);
    const auto contactIndexPtr = std::find(neighborContacts.begin(),neighborContacts.end(),uniqueIndex);
    const int  storageContactIndex = std::distance(neighborContacts.begin(),contactIndexPtr);

    storageLocation.nodeListIndex = storageNodeListIndex;
    storageLocation.nodeIndex = storageNodeIndex;
    storageLocation.contactIndex = storageContactIndex;
}

}