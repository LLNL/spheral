namespace Spheral {

//------------------------------------------------------------------------------
// set/get our spring constant
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
normalSpringConstant() const {
  return mNormalSpringConstant;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
normalSpringConstant(typename Dimension::Scalar x) {
  mNormalSpringConstant = x;
}


//------------------------------------------------------------------------------
// set/get the restitution coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
restitutionCoefficient() const {
  return mRestitutionCoefficient;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
restitutionCoefficient(typename Dimension::Scalar x) {
  mRestitutionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get this convienent damping parameter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
beta() const {
  return mBeta;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
beta(typename Dimension::Scalar x) {
  mBeta = x;
}

//------------------------------------------------------------------------------
// set/get the time step
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
timeStep() const {
  return mTimeStep;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
timeStep(typename Dimension::Scalar x) {
  mTimeStep = x;
}
}