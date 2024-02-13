namespace Spheral {

//------------------------------------------------------------------------------
// Return the set of flaw activation energies for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>
TensorDamageModel<Dimension>::
flawsForNode(const size_t index) const {
  return mFlaws(index);
}

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
TensorDamageModel<Dimension>::
youngsModulus() const {
  return mYoungsModulus;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
TensorDamageModel<Dimension>::
longitudinalSoundSpeed() const {
  return mLongitudinalSoundSpeed;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
TensorDamageModel<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
TensorDamageModel<Dimension>::
effectiveStrain() const {
  return mEffectiveStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
TensorDamageModel<Dimension>::
DdamageDt() const {
  return mDdamageDt;
}

template<typename Dimension>
inline
const typename TensorDamageModel<Dimension>::FlawStorageType&
TensorDamageModel<Dimension>::
flaws() const {
  return mFlaws;
}

template<typename Dimension>
inline
typename TensorDamageModel<Dimension>::FlawStorageType&
TensorDamageModel<Dimension>::
flaws() {
  return mFlaws;
}

template<typename Dimension>
inline
void
TensorDamageModel<Dimension>::
flaws(const typename TensorDamageModel<Dimension>::FlawStorageType& x) {
  mFlaws = x;
}

//------------------------------------------------------------------------------
// The strain update algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
TensorStrainAlgorithm
TensorDamageModel<Dimension>::
strainAlgorithm() const {
  return mStrainAlgorithm;
}

//------------------------------------------------------------------------------
// Flag to determine if damage in compression is allowed.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
TensorDamageModel<Dimension>::
damageInCompression() const {
  return mDamageInCompression;
}

template<typename Dimension>
inline
void
TensorDamageModel<Dimension>::
damageInCompression(bool x) {
  mDamageInCompression = x;
}

//------------------------------------------------------------------------------
// Threshold for the damage beyond which a node no longer votes on the time
// step.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TensorDamageModel<Dimension>::
criticalDamageThreshold() const {
  return mCriticalDamageThreshold;
}

template<typename Dimension>
inline
void
TensorDamageModel<Dimension>::
criticalDamageThreshold(double x) {
  mCriticalDamageThreshold = x;
}

}
