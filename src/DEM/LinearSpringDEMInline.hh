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
restitutionCoefficient() const {
  return mRestitutionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
restitutionCoefficient(typename Dimension::Scalar x) {
  mRestitutionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get this convienent damping parameter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
beta() const {
  return mBeta;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
beta(typename Dimension::Scalar x) {
  mBeta = x;
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