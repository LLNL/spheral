#include "SmoothingScaleBase.hh"
#include "NodeList.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace NodeSpace {

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
SmoothingScaleBase(const SmoothingScaleBase<Dimension>& rhs) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
SmoothingScaleBase<Dimension>&
SmoothingScaleBase<Dimension>::
operator=(const SmoothingScaleBase& rhs) {
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
newSmoothingScaleAndDerivative(const FieldSpace::Field<Dimension, SymTensor>& H,
                               const FieldSpace::Field<Dimension, Vector>& position,
                               const FieldSpace::Field<Dimension, Tensor>& DvDx,
                               const FieldSpace::Field<Dimension, Scalar>& zerothMoment,
                               const FieldSpace::Field<Dimension, SymTensor>& secondMoment,
                               const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                               const KernelSpace::TableKernel<Dimension>& W,
                               const Scalar hmin,
                               const Scalar hmax,
                               const Scalar hminratio,
                               const Scalar nPerh,
                               FieldSpace::Field<Dimension, SymTensor>& DHDt,
                               FieldSpace::Field<Dimension, SymTensor>& Hideal) const {
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
}

