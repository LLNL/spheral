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
currentFlaw() const {
  return mCurrentFlaw;
}

}
