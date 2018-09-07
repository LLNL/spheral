#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the total potential energy calculated in the last evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
extraEnergy() const {
  return mPotentialEnergy;
}

//------------------------------------------------------------------------------
// The mass density as a function of radius.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
massDensity(typename Dimension::Scalar r) const {
  REQUIRE(r >= 0.0);
  CHECK(mRs > 0.0);
  const Scalar r0 = r/mRs;
  const Scalar r1 = 1.0 + r0;
  return mCriticalDensity*mDeltac/(r0*r1*r1);
}

//------------------------------------------------------------------------------
// The enclosed mass as a function of radius.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
enclosedMass(typename Dimension::Scalar r) const {
  REQUIRE(r >= 0.0);
  CHECK(mRs > 0.0);
  const Scalar r0 = r/mRs;
  return 4.0/3.0*M_PI*mDeltac*mCriticalDensity*mRs*mRs*mRs*mRs*
    (1.0 + r0 - 2.0*log(1.0 + r0) - 1.0/(1.0 + r0));
}

//------------------------------------------------------------------------------
// The orbital velocity as a function of radius.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
orbitalVelocity(typename Dimension::Scalar r) const {
  REQUIRE(r >= 0.0);
  CHECK(mRs > 0.0);
  r = std::max(1.0e-10, r);
  CHECK(r > 0.0);
  const Scalar r0 = r/mRs;
  const Scalar r1 = 1.0 + r0;
  const Scalar rinverse = r/(r*r + 1.0e-20);
  const Scalar v2 = 4.0/3.0*M_PI*mDeltac*mCriticalDensity*mRs*mRs*mRs*mRs*
    mConstants.G()*((1.0/r1 + 2.0*log(r1) - 1.0)*rinverse + 1.0/(r1*r1) - 2.0/r1);
  return sqrt(v2);
}

//------------------------------------------------------------------------------
// The characteristic density.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
characteristicDensity() const {
  return mDeltac;
}

template<typename Dimension>
inline
void
NFWPotential<Dimension>::
setCharacteristicDensity(const typename Dimension::Scalar x) {
  REQUIRE(x > 0.0);
  mDeltac = x;
}

//------------------------------------------------------------------------------
// The scale radius.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
scaleRadius() const {
  return mRs;
}

template<typename Dimension>
inline
void
NFWPotential<Dimension>::
setScaleRadius(const typename Dimension::Scalar x) {
  REQUIRE(x >= 0.0);
  mRs = x;
}

//------------------------------------------------------------------------------
// The hubble constant (in units of 100 km/sec/Mpc
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
h0() const {
  return mh0;
}

template<typename Dimension>
inline
void
NFWPotential<Dimension>::
seth0(const typename Dimension::Scalar x) {
  REQUIRE(x > 0.0);
  mh0 = x;

  // Set the critical density, which is 3 H0^2 / (8 pi G).
  const Scalar H0 = 100.0 * 1.0e5 / 3.086e24 * h0() * mConstants.unitTimeSec();
  mCriticalDensity = 3.0*H0*H0/(8.0 * M_PI * mConstants.G());
}

//------------------------------------------------------------------------------
// The critical (closure) density.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
criticalDensity() const {
  return mCriticalDensity;
}

//------------------------------------------------------------------------------
// Access the origin of the profile.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
NFWPotential<Dimension>::
origin() const {
  return mOrigin;
}

template<typename Dimension>
inline
void
NFWPotential<Dimension>::
setOrigin(const typename Dimension::Vector& origin) {
  mOrigin = origin;
}

//------------------------------------------------------------------------------
// The maximum allowed fractional change in a particles potential (for 
// setting the timestep).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NFWPotential<Dimension>::
deltaPotentialFraction() const {
  return mDeltaPhiFraction;
}

template<typename Dimension>
inline
void
NFWPotential<Dimension>::
setDeltaPotentialFraction(const typename Dimension::Scalar deltaPhi) {
  REQUIRE(deltaPhi > 0.0);
  mDeltaPhiFraction = deltaPhi;
}

}
