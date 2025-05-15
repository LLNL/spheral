//---------------------------------Spheral++----------------------------------//
// SimpleMatrixMap
//
// Represents the data needed to define an explicit matrix for multiplication
// or inversion. Mainly for use in SimpleMatrixData.
//----------------------------------------------------------------------------//
#include "SimpleMatrixMap.hh"
#include "Utilities/DBC.hh"
#include "Distributed/allReduce.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
SimpleMatrixMap::
SimpleMatrixMap(int numLocalElements):
  mNumLocalElements(numLocalElements) {
  initializeIndices();
  mIndicesSet = false;
}

//------------------------------------------------------------------------------
// Get global indices
//------------------------------------------------------------------------------
void
SimpleMatrixMap::
initializeIndices() {
  // Get the cumulative sum of indices
  int sumNumIndices = distScan(mNumLocalElements, SPHERAL_OP_SUM);

  // Get the number of global elements
  mNumGlobalElements = allReduce(mNumLocalElements, SPHERAL_OP_SUM);

  // Assign the global index
  mFirstGlobalIndex = sumNumIndices - mNumLocalElements;
}

//------------------------------------------------------------------------------
// Return the calculated first global index
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
firstGlobalIndex() const {
  return mFirstGlobalIndex;
}


//------------------------------------------------------------------------------
// Return the calculated last global index
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
lastGlobalIndex() const {
  return mFirstGlobalIndex + mNumLocalElements - 1;
}

//------------------------------------------------------------------------------
// Return the local number of elements
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
numLocalElements() const {
  return mNumLocalElements;
}

//------------------------------------------------------------------------------
// Return the number of global elements
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
numGlobalElements() const {
  return mNumGlobalElements;
}

//------------------------------------------------------------------------------
// Return the number of elements per row
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
numElementsPerRow(int localIndex) const {
  CHECK(mIndicesSet);
  return (*mIndices)[localIndex].size();
}
// std::vector<int>
// SimpleMatrixMap::
// numElementsPerRow() const {
//   return mNumElementsPerRow;
// }

//------------------------------------------------------------------------------
// Return the global index, given a local index
//------------------------------------------------------------------------------
int
SimpleMatrixMap::
getGlobalIndex(int localIndex) const {
  return mFirstGlobalIndex + localIndex;
}

//------------------------------------------------------------------------------
// Check the class invariants
//------------------------------------------------------------------------------
void
SimpleMatrixMap::
checkClassInvariants() const {
}

//------------------------------------------------------------------------------
// Set the number of elements per row
//------------------------------------------------------------------------------
void
SimpleMatrixMap::
setIndices(std::vector<std::vector<int>>& indices) {
  mIndices = &indices;
  mIndicesSet = true;
}

//------------------------------------------------------------------------------
// Get nonzero columns corresponding to this row
//------------------------------------------------------------------------------
void
SimpleMatrixMap::
getColumnIndices(const int localRowIndex,
                 std::vector<int>& globalColumnIndices) const {
  CHECK(mIndicesSet);
  CHECK(mIndices->size() > size_t(localRowIndex));
  globalColumnIndices = (*mIndices)[localRowIndex];
}

} // end namespace Spheral
  
