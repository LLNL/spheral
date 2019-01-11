#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// The node extent.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
ParmetisRedistributeNodes<Dimension>::
normalizedNodeExtent() const {
  return mNormalizedNodeExtent;
}

template<typename Dimension>
inline
void
ParmetisRedistributeNodes<Dimension>::
setNormalizedNodeExtent(double extent) {
  CHECK(extent > 0.0);
  mNormalizedNodeExtent = extent;
}

}
