#include "Utilities/SpheralFunctions.hh"  // sgn
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalKernelOslo::operator==(const SphericalKernelOslo& rhs) const {
  return ((mInterp == rhs.mInterp) and
          (metamax == rhs.metamax));
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const typename SphericalKernelOslo::InterpolatorType&
SphericalKernelOslo::Winterpolator() const {
  return mInterp;
}

inline
const TableKernel<Dim<3>>&
SphericalKernelOslo::baseKernel3d() const {
  return mBaseKernel3d;
}

inline
const TableKernel<Dim<1>>&
SphericalKernelOslo::baseKernel1d() const {
  return mBaseKernel1d;
}

inline
double
SphericalKernelOslo::etamax() const {
  return metamax;
}

inline
bool
SphericalKernelOslo::useInterpolation() const {
  return mUseInterpolation;
}


inline
void
SphericalKernelOslo::useInterpolation(const bool x) {
  mUseInterpolation = x;
}

}
