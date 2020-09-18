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

//------------------------------------------------------------------------------
// Compute the time derivative of a full field of H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SmoothingScaleBase<Dimension>::
newSmoothingScaleAndDerivative(const Field<Dimension, SymTensor>& H,
                               const Field<Dimension, Vector>& position,
                               const Field<Dimension, Tensor>& DvDx,
                               const Field<Dimension, Scalar>& zerothMoment,
                               const Field<Dimension, SymTensor>& secondMoment,
                               const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const Scalar hmin,
                               const Scalar hmax,
                               const Scalar hminratio,
                               const Scalar nPerh,
                               Field<Dimension, SymTensor>& DHDt,
                               Field<Dimension, SymTensor>& Hideal) const {
  const NodeList<Dimension>& nodeList = H.nodeList();
  REQUIRE(DvDx.nodeListPtr() == &nodeList);
  REQUIRE(zerothMoment.nodeListPtr() == &nodeList);
  REQUIRE(secondMoment.nodeListPtr() == &nodeList);
  REQUIRE(DHDt.nodeListPtr() == &nodeList);
  REQUIRE(Hideal.nodeListPtr() == &nodeList);
  const unsigned nodeListi = connectivityMap.nodeListIndex(&nodeList);
  const unsigned n = nodeList.numInternalNodes();
  for (unsigned i = 0; i != n; ++i) {
    DHDt(i) = smoothingScaleDerivative(H(i),
                                       position(i),
                                       DvDx(i),
                                       hmin,
                                       hmax,
                                       hminratio,
                                       nPerh);
    Hideal(i) = newSmoothingScale(H(i), 
                                  position(i),
                                  zerothMoment(i),
                                  secondMoment(i),
                                  W,
                                  hmin,
                                  hmax,
                                  hminratio,
                                  nPerh,
                                  connectivityMap,
                                  nodeListi,
                                  i);
  }
}

}
