namespace Spheral {

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
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
const Field<Dimension, typename Dimension::SymTensor>&
TensorDamageModel<Dimension>::
newEffectiveDamage() const {
  return mNewEffectiveDamage;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Vector>&
TensorDamageModel<Dimension>::
newDamageGradient() const {
  return mNewDamageGradient;
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
// The effective damage update algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
EffectiveDamageAlgorithm
TensorDamageModel<Dimension>::
effectiveDamageAlgorithm() const {
  return mEffDamageAlgorithm;
}

//------------------------------------------------------------------------------
// Flag to determine if we compute the gradient of the damage at the start 
// of a timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
TensorDamageModel<Dimension>::
useDamageGradient() const {
  return mUseDamageGradient;
}

template<typename Dimension>
inline
void
TensorDamageModel<Dimension>::
useDamageGradient(const bool x) {
  mUseDamageGradient = x;
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
damageInCompression(const bool x) {
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
criticalDamageThreshold(const double x) {
  mCriticalDamageThreshold = x;
}

}
