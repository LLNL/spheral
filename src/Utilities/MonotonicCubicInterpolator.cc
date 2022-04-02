//---------------------------------Spheral++----------------------------------//
// MonotonicCubicInterpolator
//
// A monotonic form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#include "MonotonicCubicInterpolator.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/safeInv.hh"

#include <algorithm>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
MonotonicCubicInterpolator::MonotonicCubicInterpolator():
  mN(),
  mXmin(),
  mXmax(),
  mXstep(),
  mVals() {
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given tabluated data
//------------------------------------------------------------------------------
void
MonotonicCubicInterpolator::initialize(const double xmin,
                                       const double xmax,
                                       const std::vector<double>& yvals) {
  mN = yvals.size();
  VERIFY2(mN > 2u, "MonotonicCubicInterpolator::initialize requires at least 3 unique values to fit");
  VERIFY2(xmax > xmin, "MonotonicCubicInterpolator::initialize requires a positive domain: [" << xmin << " " << xmax << "]");

  mXmin = xmin;
  mXmax = xmax;
  mXstep = (xmax - xmin)/(mN - 1u);
  mVals.resize(2u*mN);

  // Copy the function values
  std::copy(yvals.begin(), yvals.end(), mVals.begin());

  // Estimate the gradients at our lattice points using the Catmull-Rom spline definition, and
  // then limit them to be monotone.
  this->estimateMonotoneGradients();
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
MonotonicCubicInterpolator::~MonotonicCubicInterpolator() {
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
bool
MonotonicCubicInterpolator::
operator==(const MonotonicCubicInterpolator& rhs) const {
  return ((mN == rhs.mN) and
          (mXmin == rhs.mXmin) and
          (mXmax == rhs.mXmax) and
          (mXstep == rhs.mXstep) and
          (mVals == rhs.mVals));
}

//------------------------------------------------------------------------------
// Construct monotonic limiting for the gradient values. This is the
// Fritch-Carlson method, as outlined at
// https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
//------------------------------------------------------------------------------
void
MonotonicCubicInterpolator::
estimateMonotoneGradients() {

  // Initialize gradients to zero
  std::fill(mVals.begin() + mN, mVals.end(), 0.0);

  // Compute the slope between tabulated values.
  std::vector<double> d(mN - 1u);
  const auto dxInv = 1.0/mXstep;
  for (auto k = 0u; k < mN - 1u; ++k) d[k] = (mVals[k + 1u] - mVals[k])*dxInv;

  // Initialize the tabulated gradient values using these mid-point estimates --
  // check for extrema as we go.
  mVals[mN] = d[0];
  mVals[2u*mN - 1u] = d[mN - 2u];
  std::vector<bool> mask(mN, false);
  for (auto k = 1u; k < mN - 1u; ++k) {
    if (not mask[k]) {
      if (mVals[k] == mVals[k + 1u]) {  // Leave these gradients as zero
        mask[k] = true;
        mask[k + 1u] = true;
      } else {
        mVals[mN + k] = 0.5*(d[k] + d[k + 1u]);
        // mVals[mN + k] = (d[k - 1u]*d[k] <= 0.0 ?
        //                  0.0 :
        //                  0.5*(d[k] + d[k + 1u]));
      }
    }
  }

  // // Apply limiting to the remaining unmasked gradients
  // for (auto k = 0u; k < mN - 2u; ++k) {
  //   if (not mask[k] and std::abs(d[k]) > 0.0) {
  //     const auto alpha = mVals[mN + k]*safeInv(d[k]);
  //     const auto beta = mVals[mN + k + 1u]*safeInv(d[k]);
  //     if (alpha < 0.0) {
  //       mVals[mN + k] = 0.0;
  //     } else if (beta < 0.0) {
  //       mVals[mN + k + 1u] = 0.0;
  //     }
  //     const auto tau = 3.0*safeInv(sqrt(alpha*alpha + beta*beta));
  //     mVals[mN + k] *= tau*alpha*d[k];
  //     mVals[mN + k + 1u] *= tau*beta*d[k];
  //   }
  // }
}

}
