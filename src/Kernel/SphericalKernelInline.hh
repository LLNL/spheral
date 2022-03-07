#include "Utilities/SpheralFunctions.hh"  // sgn
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalKernel::operator==(const SphericalKernel& rhs) const {
  return ((mInterp == rhs.mInterp) and
          (metamax == rhs.metamax));
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const typename SphericalKernel::InterpolatorType&
SphericalKernel::Winterpolator() const {
  return mInterp;
}

inline
const TableKernel<Dim<3>>&
SphericalKernel::baseKernel3d() const {
  return mBaseKernel3d;
}

inline
const TableKernel<Dim<1>>&
SphericalKernel::baseKernel1d() const {
  return mBaseKernel1d;
}

inline
double
SphericalKernel::etamax() const {
  return metamax;
}

inline
bool
SphericalKernel::useInterpolation() const {
  return mUseInterpolation;
}


inline
void
SphericalKernel::useInterpolation(const bool x) {
  mUseInterpolation = x;
}

}
