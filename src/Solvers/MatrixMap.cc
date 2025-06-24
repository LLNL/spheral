//---------------------------------Spheral++----------------------------------//
// MatrixMap
//
// Contains the information needed to initialize a matrix and linear solver.
// Assumes that the global indexing is sequential on one process and
// contiguous between processors. 
//----------------------------------------------------------------------------//
#include "MatrixMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
MatrixMap::
MatrixMap() {
}

int
MatrixMap::
getGlobalIndex(int localIndex) const {
  return localIndex + firstGlobalIndex();
}

std::vector<int>
MatrixMap::
numElementsPerRow() const {
  auto const numElements = this->numLocalElements();
  std::vector<int> result(numElements);
  for (auto i = 0; i < numElements; ++i) {
    result[i] = this->numElementsPerRow(i);
  }
  
  return result;
}

} // end namespace Spheral
