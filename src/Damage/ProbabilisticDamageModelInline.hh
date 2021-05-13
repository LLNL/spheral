namespace Spheral {

//------------------------------------------------------------------------------
// Access the state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
TensorStrainAlgorithm
ProbabilisticDamageModel<Dimension>::
strainAlgorithm() const {
  return mStrainAlgorithm;
}

template<typename Dimension>
inline
bool
ProbabilisticDamageModel<Dimension>::
damageInCompression() const {
  return mDamageInCompression;
}

template<typename Dimension>
inline
double
ProbabilisticDamageModel<Dimension>::
kWeibull() const {
  return mkWeibull;
}

template<typename Dimension>
inline
double
ProbabilisticDamageModel<Dimension>::
mWeibull() const {
  return mmWeibull;
}

template<typename Dimension>
inline
size_t
ProbabilisticDamageModel<Dimension>::
seed() const {
  return mSeed;
}

template<typename Dimension>
inline
size_t
ProbabilisticDamageModel<Dimension>::
minFlawsPerNode() const {
  return mMinFlawsPerNode;
}

template<typename Dimension>
inline
const Field<Dimension, unsigned>&
ProbabilisticDamageModel<Dimension>::
numFlaws() const {
  return mNumFlaws;
}

template<typename Dimension>
inline
const Field<Dimension, unsigned>&
ProbabilisticDamageModel<Dimension>::
numFlawsActivated() const {
  return mNumFlawsActivated;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
ProbabilisticDamageModel<Dimension>::
initialVolume() const {
  return mInitialVolume;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
ProbabilisticDamageModel<Dimension>::
currentFlaw() const {
  return mCurrentFlaw;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
ProbabilisticDamageModel<Dimension>::
youngsModulus() const {
  return mYoungsModulus;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
ProbabilisticDamageModel<Dimension>::
longitudinalSoundSpeed() const {
  return mLongitudinalSoundSpeed;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
ProbabilisticDamageModel<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::SymTensor>&
ProbabilisticDamageModel<Dimension>::
effectiveStrain() const {
  return mEffectiveStrain;
}

template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
ProbabilisticDamageModel<Dimension>::
DdamageDt() const {
  return mDdamageDt;
}

template<typename Dimension>
const Field<Dimension, uniform_random_01>&
ProbabilisticDamageModel<Dimension>::
randomGenerators() const {
  return mRandomGenerators;
}

}
