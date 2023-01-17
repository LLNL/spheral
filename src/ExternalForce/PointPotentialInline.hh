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
  return mTotalPotentialEnergy;
}

//------------------------------------------------------------------------------
// The speicific potential.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
specificPotential(const typename Dimension::Vector& r) const {
  return -mG*mMass/std::sqrt((mMetric*(r - mOrigin)).magnitude2() + mCoreRadius2);
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
G(const typename Dimension::Scalar x) {
  REQUIRE(x > 0.0);
  mG = x;
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
mass(const typename Dimension::Scalar x) {
  REQUIRE(x >= 0.0);
  mMass = x;
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
coreRadius(const typename Dimension::Scalar x) {
  REQUIRE(x >= 0.0);
  mCoreRadius2 = x*x;
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
origin(const typename Dimension::Vector& x) {
  mOrigin = x;
}

//------------------------------------------------------------------------------
// The metric for mapping relative positions with respect to the origin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Tensor&
PointPotential<Dimension>::
metric() const {
  return mMetric;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
metric(const typename Dimension::Tensor& x) {
  mMetric = x;
}

//------------------------------------------------------------------------------
// Timestep control
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PointPotential<Dimension>::
ftimestep() const {
  return mftimestep;
}

template<typename Dimension>
inline
void
PointPotential<Dimension>::
ftimestep(const typename Dimension::Scalar x) {
  REQUIRE(x > 0.0);
  mftimestep = x;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
PointPotential<Dimension>::
potential() const {
  return mPotential;
}

}
