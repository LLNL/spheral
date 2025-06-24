#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
inline
MatrixData::
MatrixData() {
}

//------------------------------------------------------------------------------
// Get shared_ptr for this object
//------------------------------------------------------------------------------
inline
std::shared_ptr<MatrixData>
MatrixData::
getPointer() {
  return this->shared_from_this();
}

//------------------------------------------------------------------------------
// By default, set initial guess to be RHS
//------------------------------------------------------------------------------
inline
bool
MatrixData::
initialGuessAvailable() const {
  return false;
}

//------------------------------------------------------------------------------
// Return error by default
//------------------------------------------------------------------------------
inline
void
MatrixData::
getInitialGuess(std::vector<double>& data) const {
  ASSERT2(false, "MatrixData::getInitialGuess method not set");
}

} // end namespace Spheral
