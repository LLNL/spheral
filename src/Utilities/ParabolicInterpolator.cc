//---------------------------------Spheral++----------------------------------//
// ParabolicInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#include "ParabolicInterpolator.hh"

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with tabulated data
//------------------------------------------------------------------------------
ParabolicInterpolator::ParabolicInterpolator(const double xmin,
                                             const double xmax,
                                             const std::vector<double>& yvals):
  mXmin(xmin),
  mXmax(xmax),
  mXstep((xmax - xmin)/(yvals.size() - 1)),
  mA(yvals.size()),
  mB(yvals.size()),
  mC(yvals.size()) {

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  const auto n = yvals.size();
  REQUIRE(n > 1);

  // Because we assume fixed dx step size, our A matrix of x elemements is constant and trivial!
  EMatrix A;
  A << 0.0, 0.0, 1.0,
       1.0, 1.0, 1.0,
       4.0, 2.0, 1.0;
  EMatrix Ainv = A.inverse();

  // Find the coefficient fits.
  for (auto i0 = 0u; i0 < n - 2; ++i0) {
    const auto i1 = i0 + 1;
    const auto i2 = i0 + 2;
    CHECK(i2 < n);
    EVector B, C;

    B << yvals[i0], yvals[i1], yvals[i2];
    C = Ainv*B;
    mA[i1] = C(2);
    mB[i1] = C(1);
    mC[i1] = C(0);
  }
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
ParabolicInterpolator::~ParabolicInterpolator() {
}

}
