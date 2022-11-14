namespace Spheral {

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
DEMBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
DEMBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}


//------------------------------------------------------------------------------
// set get for numerical switches
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DEMBase<Dimension>::
firstCycle() const {
  return mFirstCycle;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
firstCycle(bool x) {
  mFirstCycle = x;
}

template<typename Dimension>
inline
int
DEMBase<Dimension>::
contactRemovalFrequency() const {
  return mContactRemovalFrequency;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
contactRemovalFrequency(int x) {
  mContactRemovalFrequency = x;
}

template<typename Dimension>
inline
int
DEMBase<Dimension>::
cycle() const {
  return mCycle;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
cycle(int x) {
  mCycle = x;
}

//------------------------------------------------------------------------------
// CFL number (ratio to estimated contact duration)
//------------------------------------------------------------------------------

template<typename Dimension>
inline
typename Dimension::Scalar
DEMBase<Dimension>::
stepsPerCollision() const {
  return mStepsPerCollision;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
stepsPerCollision(typename Dimension::Scalar x) {
  mStepsPerCollision = x;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
DEMBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
DEMBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
DEMBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename DEMDimension<Dimension>::AngularVector>&
DEMBase<Dimension>::
DomegaDt() const {
  return mDomegaDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename DEMDimension<Dimension>::AngularVector>&
DEMBase<Dimension>::
omega() const {
  return mOmega;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
DEMBase<Dimension>::
uniqueIndices() const {
  return mUniqueIndices;
}


//------------------------------------------------------------------------------
// Pair things
//------------------------------------------------------------------------------


template<typename Dimension>
inline
const FieldList<Dimension, vector<int>>&
DEMBase<Dimension>::
isActiveContact() const {
  return mIsActiveContact;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<int>>&
DEMBase<Dimension>::
neighborIndices() const {
  return mNeighborIndices;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
shearDisplacement() const {
  return mShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
rollingDisplacement() const {
  return mRollingDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
torsionalDisplacement() const {
  return mTorsionalDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
DDtShearDisplacement() const {
  return mDDtShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
newShearDisplacement() const {
  return mNewShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
DDtRollingDisplacement() const {
  return mDDtRollingDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
newRollingDisplacement() const {
  return mNewRollingDisplacement;
}


template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
DDtTorsionalDisplacement() const {
  return mDDtTorsionalDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
newTorsionalDisplacement() const {
  return mNewTorsionalDisplacement;
}


template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
equilibriumOverlap() const {
  return mEquilibriumOverlap;
}

template<typename Dimension>
inline
const vector<ContactIndex>&
DEMBase<Dimension>::
contactStorageIndices() const {
  return mContactStorageIndices;
}

//------------------------------------------------------------------------------
// moment of interia specializations
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Scalar
DEMBase<Dim<1>>::
momentOfInertia(const Dim<1>::Scalar m, const Dim<1>::Scalar R) const {
  return 0.5*m*R*R;
}

template<>
inline
Dim<2>::Scalar
DEMBase<Dim<2>>::
momentOfInertia(const Dim<2>::Scalar m, const Dim<2>::Scalar R) const {
  return 0.5*m*R*R;
}

template<>
inline
Dim<3>::Scalar
DEMBase<Dim<3>>::
momentOfInertia(const Dim<3>::Scalar m, const Dim<3>::Scalar R) const {
  return 0.4*m*R*R;
}

//------------------------------------------------------------------------------
// torsion specializations
//------------------------------------------------------------------------------
template<>
inline
DEMDimension<Dim<1>>::AngularVector
DEMBase<Dim<1>>::
torsionMoment(const Dim<1>::Vector rhatij, 
              const DEMDimension<Dim<1>>::AngularVector omegai,
              const DEMDimension<Dim<1>>::AngularVector omegaj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<2>>::AngularVector
DEMBase<Dim<2>>::
torsionMoment(const Dim<2>::Vector rhatij, 
              const DEMDimension<Dim<2>>::AngularVector omegai,
              const DEMDimension<Dim<2>>::AngularVector omegaj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<3>>::AngularVector
DEMBase<Dim<3>>::
torsionMoment(const Dim<3>::Vector rhatij, 
              const DEMDimension<Dim<3>>::AngularVector omegai,
              const DEMDimension<Dim<3>>::AngularVector omegaj) const {
  return rhatij;
}


//------------------------------------------------------------------------------
// rolling Moment specializations
//------------------------------------------------------------------------------
template<>
inline
DEMDimension<Dim<1>>::AngularVector
DEMBase<Dim<1>>::
rollingMoment(const Dim<1>::Vector rhatij, 
              const Dim<1>::Vector vroti,
              const Dim<1>::Vector vrotj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<2>>::AngularVector
DEMBase<Dim<2>>::
rollingMoment(const Dim<2>::Vector rhatij, 
              const Dim<2>::Vector vroti,
              const Dim<2>::Vector vrotj) const {
  return ( DEMDimension<Dim<2>>::cross((vroti + vrotj),rhatij) > 0.0 ? 1.0 : -1.0 );
}

template<>
inline
DEMDimension<Dim<3>>::AngularVector
DEMBase<Dim<3>>::
rollingMoment(const Dim<3>::Vector rhatij, 
              const Dim<3>::Vector vroti,
              const Dim<3>::Vector vrotj) const {
  return DEMDimension<Dim<3>>::cross((vroti + vrotj),rhatij).unitVector();
}



}
