#include "Kernel/TableKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Do we require ghost-ghost or intersection connectivity?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DamageModel<Dimension>::
requireGhostConnectivity() const {
  return mDamageCouplingAlgorithm == DamageCouplingAlgorithm::ThreePointDamage;
}

template<typename Dimension>
inline
bool
DamageModel<Dimension>::
requireIntersectionConnectivity() const {
  return mComputeIntersectConnectivity;
}

//------------------------------------------------------------------------------
// Access the internal parameters of the model.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SolidNodeList<Dimension>&
DamageModel<Dimension>::
nodeList() {
  return mNodeList;
}

template<typename Dimension>
inline
const SolidNodeList<Dimension>&
DamageModel<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const TableKernel<Dimension>&
DamageModel<Dimension>::
kernel() const {
  return mW;
}

template<typename Dimension>
inline
double
DamageModel<Dimension>::
crackGrowthMultiplier() const {
  return mCrackGrowthMultiplier;
}

template<typename Dimension>
inline
DamageCouplingAlgorithm
DamageModel<Dimension>::
damageCouplingAlgorithm() const {
  return mDamageCouplingAlgorithm;
}

template<typename Dimension>
inline
const NodeCoupling&
DamageModel<Dimension>::
nodeCoupling() const {
  return *mNodeCouplingPtr;
}

}
