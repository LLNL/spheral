//---------------------------------Spheral++----------------------------------//
// SimpleMatrixData
//
// Represents the data needed to define an explicit matrix for multiplication
// or inversion.
//----------------------------------------------------------------------------//
#include "SimpleMatrixData.hh"
// #include "Distributed/Communicator.hh"
#include "SimpleMatrixMap.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
SimpleMatrixData::
SimpleMatrixData(std::shared_ptr<SimpleMatrixMap> map):
  mMap(map) {
  mIndices.resize(mMap->numLocalElements());
  mValues.resize(mMap->numLocalElements());
  mMap->setIndices(mIndices);
}

//------------------------------------------------------------------------------
// Set the data
//------------------------------------------------------------------------------
void
SimpleMatrixData::
setRowValues(int localRowIndex,
             std::vector<int> const& globalColumnIndices,
             std::vector<double> const& globalColumnValues) {
  // Check sizes
  CHECK(localRowIndex < mMap->numLocalElements());
  auto numElementsPerRow = globalColumnIndices.size();
  CONTRACT_VAR(numElementsPerRow);
  CHECK(globalColumnValues.size() == numElementsPerRow);

  // Set the indices and values
  mIndices[localRowIndex] = globalColumnIndices;
  mValues[localRowIndex] = globalColumnValues;
}

//------------------------------------------------------------------------------
// Return map, which is also a parent of this class
//------------------------------------------------------------------------------
std::shared_ptr<MatrixMap>
SimpleMatrixData::
getMap() const {
  return mMap;
}

//------------------------------------------------------------------------------
// Return the given row indices and values
//------------------------------------------------------------------------------
void
SimpleMatrixData::
getRowValues(int localRowIndex,
             std::vector<int>& globalColumnIndices,
             std::vector<double>& globalColumnValues) const {
  CHECK(localRowIndex < mMap->numLocalElements());
  globalColumnIndices = mIndices[localRowIndex];
  globalColumnValues = mValues[localRowIndex];
  CHECK(globalColumnIndices.size() == size_t(mMap->numElementsPerRow(localRowIndex)));
  CHECK(globalColumnValues.size() == size_t(mMap->numElementsPerRow(localRowIndex)));
}

//------------------------------------------------------------------------------
// Check the class invariants
//------------------------------------------------------------------------------
void
SimpleMatrixData::
SimpleMatrixData::
checkClassInvariants() const {
}


} // end namespace Spheral
  
