

//---------------------------------Spheral++----------------------------------//
// DEM type alias for the particle rotation
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContactStorageLocation_hh__
#define __Spheral_ContactStorageLocation_hh__

namespace Spheral{

  struct ContactIndex {
    ContactIndex();
    ContactIndex(const int storeNodeListIndex, 
                 const int storeNodeIndex, 
                 const int storeContactIndex,
                 const int pairNodeListIndex,
                 const int pairNodeIndex);

      int storeNodeList; 
      int storeNode; 
      int storeContact;

      int pairNodeList;
      int pairNode;
  };


#else

// Forward declaration.
namespace Spheral {
  struct ContactIndex;
}

#endif
}