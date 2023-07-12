#include "Utilities/SpheralFunctions.hh"  // sgn
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalRadialKernel::operator==(const SphericalRadialKernel& rhs) const {
  return ((mAInterp == rhs.mAInterp) and
          (mGradAInvInterp == rhs.mGradAInvInterp) and
          (mBaseKernel == rhs.mBaseKernel) and
          (mNumIntegral == rhs.mNumIntegral) and
          (mUseInterpolation == rhs.mUseInterpolation));
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const typename SphericalRadialKernel::InterpolatorType&
SphericalRadialKernel::Ainterpolator() const {
  return mAInterp;
}

inline
const typename SphericalRadialKernel::InterpolatorType&
SphericalRadialKernel::gradAInvInterpolator() const {
  return mGradAInvInterp;
}

inline
const TableKernel<Dim<1>>&
SphericalRadialKernel::baseKernel1d() const {
  return mBaseKernel;
}

inline
double
SphericalRadialKernel::etamax() const {
  return metamax;
}

inline
bool
SphericalRadialKernel::useInterpolation() const {
  return mUseInterpolation;
}


inline
void
SphericalRadialKernel::useInterpolation(const bool x) {
  mUseInterpolation = x;
}

}
