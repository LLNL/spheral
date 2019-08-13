//------------------------------------------------------------------------------
// A struct for encoding cell face information returned by computeVoronoiVolume
//------------------------------------------------------------------------------
#ifndef __Spheral_CellFaceFlags__
#define __Spheral_CellFaceFlags__

#include <vector>
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/packElement.hh"

namespace Spheral {

struct CellFaceFlags {
  int cellFace;  // The index of the face in the cell.facets array
  int nodeListj; // The NodeList of the opposite node across the face
  int j;         // The index of the opposite node across the face
};

//------------------------------------------------------------------------------
// Communication/serialization traits
//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<CellFaceFlags> {
  typedef int ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const CellFaceFlags& x) { return 3; }
  static CellFaceFlags zero() { return CellFaceFlags({-1, -1, -1}); }
};

template<>
void packElement<CellFaceFlags>(const CellFaceFlags& value,
                                std::vector<char>& buffer) {
  packElement(value.cellFace, buffer);
  packElement(value.nodeListj, buffer);
  packElement(value.j, buffer);
}

template<>
void unpackElement<CellFaceFlags>(CellFaceFlags& value,
                                  std::vector<char>::const_iterator& itr,
                                  const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.cellFace, itr, endPackedVector);
  unpackElement(value.nodeListj, itr, endPackedVector);
  unpackElement(value.j, itr, endPackedVector);
  ENSURE(itr <= endPackedVector);
}

}

#endif
