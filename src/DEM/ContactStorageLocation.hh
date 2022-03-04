

//---------------------------------Spheral++----------------------------------//
// DEM type alias for the particle rotation
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContactStorageLocation_hh__
#define __Spheral_ContactStorageLocation_hh__

namespace Spheral{

  struct ContactIndex {
    ContactIndex(const int nodeListIndex, 
                 const int nodeIndex, 
                 const int contactIndex);

    int nodeListIndex; 
    int nodeIndex; 
    int contactIndex;
  }


  struct ContactStorageLocation{

    ContactStorageLocation(const int nodeListi,
                          const int i,
                          const int uniqueIndexi,
                          const int numInternalElementsi,
                          const int nodeListj,
                          const int j,
                          const int uniqueIndexj,
                          const int numInternalElementsj);
    ContactIndex storageLocation;
    ContactIndex pairLocation;               
  }
}


#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> struct ContactIndex;
  template<typename Dimension> struct ContactStorageLocation;
}

#endif
}