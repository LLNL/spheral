#include "SmoothingScaleBase.hh"
#include "NodeList.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>::
SmoothingScaleBase() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>::
SmoothingScaleBase(const SmoothingScaleBase<Dimension>& ) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>&
SmoothingScaleBase<Dimension>::
operator=(const SmoothingScaleBase&) {
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>::
~SmoothingScaleBase<Dimension>() {
}

}
