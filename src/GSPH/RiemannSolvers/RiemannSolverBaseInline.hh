namespace Spheral {

//------------------------------------------------------------------------------
// Get slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
LimiterBase<Dimension>&
RiemannSolverBase<Dimension>::
limiter() const {
  return mSlopeLimiter;
}

//------------------------------------------------------------------------------
// Get wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
inline
WaveSpeedBase<Dimension>&
RiemannSolverBase<Dimension>::
waveSpeed() const {
  return mWaveSpeed;
}

//------------------------------------------------------------------------------
// set/get linear reconstruction switch
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RiemannSolverBase<Dimension>::
linearReconstruction() const {
  return mLinearReconstruction;
}

template<typename Dimension>
inline
void
RiemannSolverBase<Dimension>::
linearReconstruction(bool x) {
  mLinearReconstruction=x;
}


//------------------------------------------------------------------------------
// field getters
//------------------------------------------------------------------------------


// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Vector>&
// RiemannSolverBase<Dimension>::
// DpDx() {
//   return mDpDx;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Tensor>&
// RiemannSolverBase<Dimension>::
// DvDx() {
//   return mDvDx;
// }


// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// RiemannSolverBase<Dimension>::
// DpDx() const {
//   return mDpDx;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Tensor>&
// RiemannSolverBase<Dimension>::
// DvDx() const {
//   return mDvDx;
// }


}