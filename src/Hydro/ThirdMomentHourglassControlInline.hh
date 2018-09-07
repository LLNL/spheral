#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the maximum multiplicative factor for the acceleration.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
ThirdMomentHourglassControl<Dimension>::
maxAccelerationFactor() const {
  return mMaxAccelerationFactor;
}

template<typename Dimension>
inline
void
ThirdMomentHourglassControl<Dimension>::
maxAccelerationFactor(const double x) {
  VERIFY(x >= 0.0);
  mMaxAccelerationFactor = x;
}

//------------------------------------------------------------------------------
// Access the multiplier.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
ThirdMomentHourglassControl<Dimension>::
multiplier() const {
  return mMultiplier;
}

template<typename Dimension>
inline
void
ThirdMomentHourglassControl<Dimension>::
multiplier(const double x) {
  VERIFY(x >= 0.0);
  mMultiplier = x;
}

//------------------------------------------------------------------------------
// The last computed third moment of the node distribution.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::ThirdRankTensor>&
ThirdMomentHourglassControl<Dimension>::
thirdMoment() const {
  return mThirdMoment;
}

}
