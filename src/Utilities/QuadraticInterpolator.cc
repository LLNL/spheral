//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#include "QuadraticInterpolator.hh"
#include "umpire/Umpire.hpp"
#include <algorithm>

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor with sampled values
//------------------------------------------------------------------------------
QIHandler::QIHandler(double xmin,
                     double xmax,
                     const std::vector<double>& yvals) {
  initialize(xmin, xmax, yvals);
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
void
QIHandler::initialize(double xmin,
                      double xmax,
                      const std::vector<double>& yvals) {
  const auto n = yvals.size();
  VERIFY2(n > 2, "QIHandler::initialize requires at least 3 unique values to fit");
  VERIFY2(n % 2 == 1, "QIHandler::initialize requires an odd number of tabulated values");
  VERIFY2(xmax > xmin, "QIHandler::initialize requires a positive domain: [" << xmin << " " << xmax << "]");

  double N1 = (n - 1u)/2u - 1u;  // Maximum index into arrays
  double xstep = (xmax - xmin)/(N1 + 1u);
  size_t N = 3*(N1 + 1u);
  auto& rm = umpire::ResourceManager::getInstance();
  auto host_allocator = rm.getAllocator("HOST");
  double* vals = static_cast<double*>(host_allocator.allocate(N*sizeof(double)));

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  // We use simple least squares fitting for each 3-point interval, giving us an
  // exact parabolic fit for those three points.
  // Find the coefficient fits.
  double x0, x1, x2;
  EMatrix A;
  EVector X, B;
  for (auto i0 = 0u; i0 <= N1; ++i0) {
    x0 = xmin + i0*xstep;
    x1 = x0 + 0.5*xstep;
    x2 = x0 + xstep;
    A << 1.0, x0, x0*x0,
         1.0, x1, x1*x1,
         1.0, x2, x2*x2;
    B << yvals[2u*i0], yvals[2u*i0 + 1u], yvals[2u*i0 + 2u];
    X = A.inverse()*B;
    vals[3*i0     ] = X(0);
    vals[3*i0 + 1u] = X(1);
    vals[3*i0 + 2u] = X(2);
  }
#ifdef SPHERAL_GPU_ENABLED
  auto allocator = rm.getAllocator("DEVICE");
  mDeviceCoeffs = static_cast<double*>(allocator.allocate(N*sizeof(double)));
  rm.copy(mDeviceCoeffs, vals, N*sizeof(double));
#endif
  mHostCoeffs = vals;
  mN1 = N1;
  mXmin = xmin;
  mXmax = xmax;
  mXstep = xstep;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
QIHandler::~QIHandler() {
  auto& rm = umpire::ResourceManager::getInstance();
#ifdef SPHERAL_GPU_ENABLED
  auto allocator = rm.getAllocator("DEVICE");
  allocator.deallocate(mDeviceCoeffs);
#endif
  auto host_allocator = rm.getAllocator("HOST");
  host_allocator.deallocate(mHostCoeffs);
}

}
