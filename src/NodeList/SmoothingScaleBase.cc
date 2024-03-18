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
                               const Field<Dimension, Vector>& firstMoment,
                               const Field<Dimension, SymTensor>& secondMomentEta,
                               const Field<Dimension, SymTensor>& secondMomentLab,
                               const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const Scalar hmin,
                               const Scalar hmax,
                               const Scalar hminratio,
                               const Scalar nPerh,
                               Field<Dimension, SymTensor>& DHDt,
                               Field<Dimension, SymTensor>& Hideal) const {
  const auto& nodeList = H.nodeList();
  REQUIRE(DvDx.nodeListPtr() == &nodeList);
  REQUIRE(zerothMoment.nodeListPtr() == &nodeList);
  REQUIRE(firstMoment.nodeListPtr() == &nodeList);
  REQUIRE(secondMomentEta.nodeListPtr() == &nodeList);
  REQUIRE(secondMomentLab.nodeListPtr() == &nodeList);
  REQUIRE(DHDt.nodeListPtr() == &nodeList);
  REQUIRE(Hideal.nodeListPtr() == &nodeList);
  const auto nodeListi = connectivityMap.nodeListIndex(&nodeList);
  const auto n = nodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
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
                                  firstMoment(i),
                                  secondMomentEta(i),
                                  secondMomentLab(i),
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
