#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// The effective critical number of nodes per smoothing scale, below which we
// assume all flaws are active on a node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
DamageModel<Dimension>::
criticalNodesPerSmoothingScale() const {
  return mCriticalNodesPerSmoothingScale;
}

template<typename Dimension>
inline
void
DamageModel<Dimension>::
criticalNodesPerSmoothingScale(const double x) {
  VERIFY(x >= 0.0);
  mCriticalNodesPerSmoothingScale = x;
}

//------------------------------------------------------------------------------
// Return the set of flaw activation energies for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>
DamageModel<Dimension>::
flawsForNode(const size_t index) const {
  REQUIRE(index < mNodeList.numInternalNodes());
  REQUIRE(mFlaws.nodeListPtr() == &mNodeList);
  if (mEffectiveFlawAlgorithm == EffectiveFlawAlgorithm::FullSpectrumFlaws) {
    return mFlaws(index);
  } else {
    return std::vector<double>(1, mEffectiveFlaws(index));
  }
}

//------------------------------------------------------------------------------
// Access the internal parameters of the model.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SolidMaterial::SolidNodeList<Dimension>&
DamageModel<Dimension>::
nodeList() {
  return mNodeList;
}

template<typename Dimension>
inline
const SolidMaterial::SolidNodeList<Dimension>&
DamageModel<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const KernelSpace::TableKernel<Dimension>&
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
EffectiveFlawAlgorithm
DamageModel<Dimension>::
effectiveFlawAlgorithm() const {
  return mEffectiveFlawAlgorithm;
}

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
DamageModel<Dimension>::
youngsModulus() const {
  return mYoungsModulus;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
DamageModel<Dimension>::
longitudinalSoundSpeed() const {
  return mLongitudinalSoundSpeed;
}

//------------------------------------------------------------------------------
// Access the flaw activation strains.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename DamageModel<Dimension>::FlawStorageType&
DamageModel<Dimension>::
flaws() const {
  return mFlaws;
}

template<typename Dimension>
inline
typename DamageModel<Dimension>::FlawStorageType&
DamageModel<Dimension>::
flaws() {
  return mFlaws;
}

//------------------------------------------------------------------------------
// Access the computed effective flaw activation strain per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
DamageModel<Dimension>::
effectiveFlaws() const {
  return mEffectiveFlaws;
}

}
}
