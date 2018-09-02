namespace Spheral {

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
TotalHydro<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
TotalHydro<Dimension>::
HEvolution(const HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Access the minimum allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
TotalHydro<Dimension>::hmin() const {
  CHECK(mhmin > 0.0);
  return mhmin;
}

template<typename Dimension>
inline
void
TotalHydro<Dimension>::
hmin(const typename Dimension::Scalar val) {
  CHECK(val > 0.0);
  mhmin = val;
}

//------------------------------------------------------------------------------
// Access the maximum allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
TotalHydro<Dimension>::hmax() const {
  CHECK(mhmax > 0.0);
  return mhmax;
}

template<typename Dimension>
inline
void
TotalHydro<Dimension>::
hmax(const typename Dimension::Scalar val) {
  CHECK(val > 0.0);
  mhmax = val;
}

//------------------------------------------------------------------------------
// Access the minimum allowed ratio of the smoothing scales in the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
TotalHydro<Dimension>::hratiomin() const {
  CHECK(mhratiomin >= 0.0);
  return mhratiomin;
}

template<typename Dimension>
inline
void
TotalHydro<Dimension>::
hratiomin(typename Dimension::Scalar val) {
  CHECK(val >= 0.0);
  mhratiomin = val;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
TotalHydro<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
TotalHydro<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
positionWeight() const {
  return mPositionWeight;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
totalEnergy() const {
  return mTotalEnergy;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
TotalHydro<Dimension>::
linearMomentum() const {
  return mLinearMomentum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
DVDt() const {
  return mDVDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
TotalHydro<Dimension>::
DEDt() const {
  return mDEDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
TotalHydro<Dimension>::
DpmomDt() const {
  return mDpmomDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
TotalHydro<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

}
