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



//------------------------------------------------------------------------------
// Pair things
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
DEMBase<Dimension>::
uniqueIndices() const {
  return mUniqueIndices;
}

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


}
