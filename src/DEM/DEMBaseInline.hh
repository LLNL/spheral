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
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
DEMBase<Dimension>::
kernel() const {
  return mKernel;
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
const FieldList<Dimension, typename Dimension::Vector>&
DEMBase<Dimension>::
DomegaDt() const {
  return mDomegaDt;
}

//------------------------------------------------------------------------------
// Access the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<ContactModelBase<Dimension>*>&
DEMBase<Dimension>::contactModels() const {
  return mContactModels;
}

//------------------------------------------------------------------------------
// Provide iterators over over the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<ContactModelBase<Dimension>*>::iterator
DEMBase<Dimension>::contactModelsBegin() {
  return mContactModels.begin();
}

template<typename Dimension>
inline
typename std::vector<ContactModelBase<Dimension>*>::iterator
DEMBase<Dimension>::contactModelsEnd() {
  return mContactModels.end();
}

//------------------------------------------------------------------------------
// Provide const iterators over over the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<ContactModelBase<Dimension>*>::const_iterator
DEMBase<Dimension>::contactModelsBegin() const {
  return mContactModels.begin();
}

template<typename Dimension>
inline
typename std::vector<ContactModelBase<Dimension>*>::const_iterator
DEMBase<Dimension>::contactModelsEnd() const {
  return mContactModels.end();
}

}
