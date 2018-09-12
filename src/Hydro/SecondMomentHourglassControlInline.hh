#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the maximum multiplicative factor for the acceleration.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SecondMomentHourglassControl<Dimension>::
maxAccelerationFactor() const {
  return mMaxAccelerationFactor;
}

template<typename Dimension>
inline
void
SecondMomentHourglassControl<Dimension>::
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
SecondMomentHourglassControl<Dimension>::
multiplier() const {
  return mMultiplier;
}

template<typename Dimension>
inline
void
SecondMomentHourglassControl<Dimension>::
multiplier(const double x) {
  VERIFY(x >= 0.0);
  mMultiplier = x;
}

//------------------------------------------------------------------------------
// Acceleration diagnostic.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SecondMomentHourglassControl<Dimension>::
acceleration() const {
  return mAcceleration;
}

}
