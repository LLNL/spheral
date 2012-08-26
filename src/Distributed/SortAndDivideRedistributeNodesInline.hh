#include "DBC.hh"
#include "Infrastructure/SpheralFunctions.hh"

namespace Spheral {
namespace PartitionSpace {

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
SortAndDivideRedistributeNodes<Dimension>::Hextent(const double val) {
  REQUIRE(distinctlyGreaterThan(val, 0.0));
  mHextent = val;
}

}
}
