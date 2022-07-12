namespace Spheral {

//------------------------------------------------------------------------------
// set/get our spring constant
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalSpringConstant() const {
  return mNormalSpringConstant;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalSpringConstant(typename Dimension::Scalar x) {
  mNormalSpringConstant = x;
}


//------------------------------------------------------------------------------
// set/get the restitution coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalRestitutionCoefficient() const {
  return mNormalRestitutionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalRestitutionCoefficient(typename Dimension::Scalar x) {
  mNormalRestitutionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get our spring constant
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialSpringConstant() const {
  return mTangentialSpringConstant;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialSpringConstant(typename Dimension::Scalar x) {
  mTangentialSpringConstant = x;
}


//------------------------------------------------------------------------------
// set/get the restitution coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialRestitutionCoefficient() const {
  return mTangentialRestitutionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialRestitutionCoefficient(typename Dimension::Scalar x) {
  mTangentialRestitutionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get the dynamic friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
dynamicFrictionCoefficient() const {
  return mDynamicFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
dynamicFrictionCoefficient(typename Dimension::Scalar x) {
  mDynamicFrictionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get the static friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
staticFrictionCoefficient() const {
  return mStaticFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
staticFrictionCoefficient(typename Dimension::Scalar x) {
  mStaticFrictionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get the rolling friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
rollingFrictionCoefficient() const {
  return mRollingFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
rollingFrictionCoefficient(typename Dimension::Scalar x) {
  mRollingFrictionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get the torsional friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
torsionalFrictionCoefficient() const {
  return mTorsionalFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
torsionalFrictionCoefficient(typename Dimension::Scalar x) {
  mTorsionalFrictionCoefficient = x;
}

//------------------------------------------------------------------------------
// cohesive coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
cohesiveCoefficient() const {
  return mCohesiveCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
cohesiveCoefficient(typename Dimension::Scalar x) {
  mCohesiveCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get the shape factor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
shapeFactor() const {
  return mShapeFactor;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
shapeFactor(typename Dimension::Scalar x) {
  mShapeFactor = x;
}

//------------------------------------------------------------------------------
// set/get this convienent damping parameter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalBeta() const {
  return mNormalBeta;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalBeta(typename Dimension::Scalar x) {
  mNormalBeta = x;
}

template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialBeta() const {
  return mTangentialBeta;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialBeta(typename Dimension::Scalar x) {
  mTangentialBeta = x;
}

//------------------------------------------------------------------------------
// set/get the time step
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
timeStep() const {
  return mTimeStep;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
timeStep(typename Dimension::Scalar x) {
  mTimeStep = x;
}


}