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
    ContactIndex(const int storeNodeListIndex, 
                 const int storeNodeIndex, 
                 const int storeContactIndex,
                 const int solidBoundaryIndex);

      int storeNodeList; 
      int storeNode; 
      int storeContact;

      int pairNodeList;
      int pairNode;

      int solidBoundary;
  };
}

#else

// Forward declaration.
namespace Spheral {
  struct ContactIndex;
}

#endif
