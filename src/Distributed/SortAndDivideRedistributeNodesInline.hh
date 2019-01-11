#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the Hextent variable.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SortAndDivideRedistributeNodes<Dimension>::Hextent() const {
  return mHextent;
}

template<typename Dimension>
inline
void
SortAndDivideRedistributeNodes<Dimension>::Hextent(double val) {
  REQUIRE(distinctlyGreaterThan(val, 0.0));
  mHextent = val;
}

}
