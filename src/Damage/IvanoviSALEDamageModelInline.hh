namespace Spheral {

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
IvanoviSALEDamageModel<Dimension>::
minPlasticFailure() const {
  return mEpsPfb;
}

template<typename Dimension>
inline
double
IvanoviSALEDamageModel<Dimension>::
plasticFailurePressureSlope() const {
  return mB;
}

template<typename Dimension>
inline
double
IvanoviSALEDamageModel<Dimension>::
plasticFailurePressureOffset() const {
  return mPc;
}

template<typename Dimension>
inline
double
IvanoviSALEDamageModel<Dimension>::
tensileFailureStress() const {
  return mTensileFailureStress;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
IvanoviSALEDamageModel<Dimension>::
youngsModulus() const {
  return mYoungsModulus;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
IvanoviSALEDamageModel<Dimension>::
longitudinalSoundSpeed() const {
  return mLongitudinalSoundSpeed;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
IvanoviSALEDamageModel<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
IvanoviSALEDamageModel<Dimension>::
effectiveStrain() const {
  return mEffectiveStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
IvanoviSALEDamageModel<Dimension>::
DdamageDt() const {
  return mDdamageDt;
}

template<typename Dimension>
const Field<Dimension, int>&
IvanoviSALEDamageModel<Dimension>::
mask() const {
  return mMask;
}

template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
mask(const Field<Dimension, int>& val) {
  mMask = val;
}

template<typename Dimension>
double
IvanoviSALEDamageModel<Dimension>::
criticalDamageThreshold() const {
  return mCriticalDamageThreshold;
}

template<typename Dimension>
void
IvanoviSALEDamageModel<Dimension>::
criticalDamageThreshold(const double val) {
  mCriticalDamageThreshold = val;
}

}
