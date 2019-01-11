#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the total potential energy calculated in the last evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
extraEnergy() const {
  return mPotentialEnergy;
}

//------------------------------------------------------------------------------
// The speicific potential.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
specificPotential(const typename Dimension::Vector& r) const {
  return -G()*mass()/sqrt((r - origin()).magnitude2() + mCoreRadius2);
}

//------------------------------------------------------------------------------
// Access the gravitational constant being used.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
G() const {
  return mG;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
setG(const typename Dimension::Scalar G) {
  REQUIRE(G > 0.0);
  mG = G;
}

//------------------------------------------------------------------------------
// Access the point mass being used.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
mass() const {
  return mMass;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
setMass(const typename Dimension::Scalar m) {
  REQUIRE(m >= 0.0);
  mMass = m;
}

//------------------------------------------------------------------------------
// Access the core radius being used.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
coreRadius() const {
  return sqrt(mCoreRadius2);
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
setCoreRadius(const typename Dimension::Scalar rc) {
  REQUIRE(rc >= 0.0);
  mCoreRadius2 = rc*rc;
}

//------------------------------------------------------------------------------
// Access the position of the point mass.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
PointPotential<Dimension>::
origin() const {
  return mOrigin;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
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
PointPotential<Dimension>::
deltaPotentialFraction() const {
  return mDeltaPhiFraction;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
setDeltaPotentialFraction(const typename Dimension::Scalar deltaPhi) {
  REQUIRE(deltaPhi > 0.0);
  mDeltaPhiFraction = deltaPhi;
}

}
