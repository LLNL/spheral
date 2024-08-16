//---------------------------------Spheral++----------------------------------//
// CubicHermiteInterpolator
//
// An (optionally monotonic) form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#include "CubicHermiteInterpolator.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/safeInv.hh"

#include <algorithm>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct from explicit table of values
//------------------------------------------------------------------------------
CubicHermiteInterpolator::
CubicHermiteInterpolator(const double xmin,
                         const double xmax,
                         const std::vector<double>& yvals):
  mN(),
  mXmin(),
  mXmax(),
  mXstep(),
  mVals() {
  this->initialize(xmin, xmax, yvals);
}

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
CubicHermiteInterpolator::CubicHermiteInterpolator():
  mN(),
  mXmin(),
  mXmax(),
  mXstep(),
  mVals() {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
CubicHermiteInterpolator::CubicHermiteInterpolator(const CubicHermiteInterpolator& rhs):
  mN(rhs.mN),
  mXmin(rhs.mXmin),
  mXmax(rhs.mXmax),
  mXstep(rhs.mXstep),
  mVals(rhs.mVals) {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
CubicHermiteInterpolator&
CubicHermiteInterpolator::operator=(const CubicHermiteInterpolator& rhs) {
  if (this != &rhs) {
    mN = rhs.mN;
    mXmin = rhs.mXmin;
    mXmax = rhs.mXmax;
    mXstep = rhs.mXstep;
    mVals = rhs.mVals;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given tabluated data
//------------------------------------------------------------------------------
void
CubicHermiteInterpolator::initialize(const double xmin,
                                     const double xmax,
                                     const std::vector<double>& yvals) {
  mN = yvals.size();
  VERIFY2(mN > 2u, "CubicHermiteInterpolator::initialize requires at least 3 unique values to fit");
  VERIFY2(xmax > xmin, "CubicHermiteInterpolator::initialize requires a positive domain: [" << xmin << " " << xmax << "]");

  mXmin = xmin;
  mXmax = xmax;
  mXstep = (xmax - xmin)/(mN - 1u);
  mVals.resize(2u*mN);

  // Copy the function values
  std::copy(yvals.begin(), yvals.end(), mVals.begin());

  // Estimate the gradients at our lattice points
  this->initializeGradientKnots();
  // const auto dxInv = 1.0/mXstep;
  // for (auto i = 1u; i < mN - 1u; ++i) {
  //   mVals[mN + i] = 0.5*(mVals[i + 1u] - mVals[i - 1u])*dxInv;
  // }
  // mVals[mN] = (mVals[1] - mVals[0])*dxInv;
  // mVals[2u*mN - 1u] = (mVals[mN - 1u] - mVals[mN - 2u])*dxInv;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
CubicHermiteInterpolator::~CubicHermiteInterpolator() {
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
bool
CubicHermiteInterpolator::
operator==(const CubicHermiteInterpolator& rhs) const {
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
CubicHermiteInterpolator::
makeMonotonic() {

  // Compute the slope between tabulated values.
  std::vector<double> cgrad(mN - 1u);
  const auto dxInv = 1.0/mXstep;
  for (auto k = 0u; k < mN - 1u; ++k) cgrad[k] = (mVals[k + 1u] - mVals[k])*dxInv;

  // Helper function for indexing the gradent values
  auto gradVal = [&](const size_t i) -> double& { return mVals[mN + i]; };

  // Check for extrema (places where the centered gradient changes sign), and zero
  // out the gradient there
  for (auto k = 1u; k < mN - 1u; ++k) {
    if (cgrad[k - 1u]*cgrad[k] <= 0.0) gradVal(k) = 0.0;
    if (cgrad[k] == 0.0) {
      gradVal(k) = 0.0;
      gradVal(k + 1u) = 0.0;
    }
  }

  // Iterative limit the slopes until we stop changing values
  auto done = false;
  while (not done) {
    done = true;
    for (auto k = 0u; k < mN - 1u; ++k) {
      auto alpha = gradVal(k)/cgrad[k];
      auto beta = gradVal(k + 1u)/cgrad[k];
      if (alpha < 0.0) {
        gradVal(k) = 0.0;
        alpha = 0.0;
        done = false;
      }
      if (beta < 0.0) {
        gradVal(k + 1u) = 0.0;
        beta = 0.0;
        done = false;
      }
      const auto tau = 3.0/sqrt(alpha*alpha + beta*beta);
      if (tau < 1.0) {
        gradVal(k) = 0.99*tau*alpha*cgrad[k];
        gradVal(k + 1u) = 0.99*tau*beta*cgrad[k];
        done = false;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Construct the gradient values at the knot points using the algorithm
// described in 
// Note: the function values must have already been set in mVals[0,n-1]!
//------------------------------------------------------------------------------
void
CubicHermiteInterpolator::
initializeGradientKnots() {

  // Set up to solve using Eigen's linear solvers (Ax=b)
  Eigen::MatrixXf A = Eigen::MatrixXf::Zero(mN, mN);
  Eigen::VectorXf b(mN);
  A(0, 0)         =  4.0;      // First row
  A(0, 1)         = -1.0;
  A(mN-1u, mN-2u) = -1.0;      // Last row
  A(mN-1u, mN-1u) =  4.0;
  b(0)            =  3.0*(mVals[1u] - mVals[0u])/mXstep;
  b(mN-1u)        =  3.0*(mVals[mN-1u] - mVals[mN-2u])/mXstep;
  for (auto k = 1u; k < mN-1u; ++k) {      // rows
    A(k, k-1u) = -0.5;
    A(k, k)    = 4.0;
    A(k, k+1u) = -0.5;
    b(k) = 1.5*(mVals[k+1u] - mVals[k-1u])/mXstep;
  }

  // Solve for the gradient values
  const Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);
  for (auto k = 0u; k < mN; ++k) mVals[mN + k] = x(k);
}

}
