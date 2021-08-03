namespace Spheral {

//------------------------------------------------------------------------------
// Set/Get slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SlopeLimiterBase<Dimension>&
RiemannSolverBase<Dimension>::
slopeLimiter() const {
  return mSlopeLimiter;
}

template<typename Dimension>
inline
void
RiemannSolverBase<Dimension>::
slopeLimiter(SlopeLimiterBase<Dimension>& slopeLimiter) {
  mSlopeLimiter=slopeLimiter;
}

//------------------------------------------------------------------------------
// Set/Get wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const WaveSpeedBase<Dimension>&
RiemannSolverBase<Dimension>::
waveSpeed() const {
  return mWaveSpeed;
}

template<typename Dimension>
inline
void
RiemannSolverBase<Dimension>::
waveSpeed(WaveSpeedBase<Dimension>& waveSpeed) {
  mWaveSpeed=waveSpeed;
}

//------------------------------------------------------------------------------
// field getters
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
RiemannSolverBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
RiemannSolverBase<Dimension>::
DpDx() const {
  return mDpDx;
}


}