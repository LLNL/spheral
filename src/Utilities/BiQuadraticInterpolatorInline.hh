#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
BiQuadraticInterpolator::BiQuadraticInterpolator():
  mnx1(),
  mny1(),
  mxlog(),
  mylog(),
  mxmin(),
  mxmax(),
  mxstep(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Construct with tabulated data
//------------------------------------------------------------------------------
template<typename Func>
BiQuadraticInterpolator::BiQuadraticInterpolator(const Vector& xmin,
                                                 const Vector& xmax,
                                                 const size_t nx,
                                                 const size_t ny,
                                                 const bool logxspace,
                                                 const bool logyspace,
                                                 const Func& F):
  mnx1(),
  mny1(),
  mxlog(),
  mylog(),
  mxmin(),
  mxmax(),
  mxstep(),
  mcoeffs() {
  this->initialize(xmin, xmax, nx, ny, logxspace, logyspace, F);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
BiQuadraticInterpolator::~BiQuadraticInterpolator() {
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
template<typename Func>
void
BiQuadraticInterpolator::initialize(const Vector& xmin,
                                    const Vector& xmax,
                                    const size_t nx,
                                    const size_t ny,
                                    const bool logxspace,
                                    const bool logyspace,
                                    const Func& F) {

  typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 6, 1> EVector;

  // Size stuff up.
  REQUIRE(nx > 1);
  REQUIRE(ny > 1);
  mnx1 = nx - 1u;
  mny1 = ny - 1u;
  mxlog = logxspace;
  mylog = logyspace;
  mcoeffs.resize(6*mnx1*mny1);

  // Figure out the sampling steps.
  mxmin = xmin;
  mxmax = xmax;
  if (logxspace) {
    mxmin[0] = log(mxmin[0]);
    mxmax[0] = log(mxmax[0]);
  }
  if (logyspace) {
    mxmin[1] = log(mxmin[1]);
    mxmax[1] = log(mxmax[1]);
    mxstep[1] = (log(xmax[1]) - log(xmin[1]))/mny1;
  }
  mxstep = {(xmax[0] - xmin[0])/mnx1,
            (xmax[1] - xmin[1])/mny1};

  // Fit the coefficients
  Vector x00, x01, x10, x11, xmid1, xmid2;
  EMatrix A;
  EVector B, C;
  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x00[0] = xmin[0] + i    *mxstep[0];
      x11[0] = xmin[0] + (i + 1)*mxstep[0];
      xmid1[0] = xmin[0] + (i + 0.3333333333)*mxstep[0];
      xmid2[0] = xmin[0] + (i + 0.6666666666)*mxstep[0];
      if (logxspace) {
        x00[0] = exp(x00[0]);
        x11[0] = exp(x11[0]);
        xmid1[0] = exp(xmid1[0]);
        xmid2[0] = exp(xmid2[0]);
      }
      x00[1] = xmin[1] + j    *mxstep[1];
      x11[1] = xmin[1] + (j + 1)*mxstep[1];
      xmid1[1] = xmin[1] + (j + 0.3333333333)*mxstep[1];
      xmid2[1] = xmin[1] + (j + 0.6666666666)*mxstep[1];
      if (logyspace) {
        x00[1] = exp(x00[1]);
        x11[1] = exp(x11[1]);
        xmid1[1] = exp(xmid1[1]);
        xmid2[1] = exp(xmid2[1]);
      }
      x01 = {x00[0], x11[1]};
      x10 = {x11[0], x00[1]};
      A << 1.0, x00[0], x00[1], x00[0]*x00[1], x00[0]*x00[0], x00[1]*x00[1],
           1.0, x01[0], x01[1], x01[0]*x01[1], x01[0]*x01[0], x01[1]*x01[1],
           1.0, x10[0], x10[1], x10[0]*x10[1], x10[0]*x10[0], x10[1]*x10[1],
           1.0, x11[0], x11[1], x11[0]*x11[1], x11[0]*x11[0], x11[1]*x11[1],
           1.0, xmid1[0], xmid1[1], xmid1[0]*xmid1[1], xmid1[0]*xmid1[0], xmid1[1]*xmid1[1],
           1.0, xmid2[0], xmid2[1], xmid2[0]*xmid2[1], xmid2[0]*xmid2[0], xmid2[1]*xmid2[1];
      B << F(x00), F(x01), F(x10), F(x11), F(xmid1), F(xmid2);
      C = A.inverse()*B;
      std::cerr << "------------------------------------------------------------------------------\n"
                << "x00: " << x00 << "\n"
                << "x10: " << x10 << "\n"
                << "x01: " << x01 << "\n"
                << "x11: " << x11 << "\n"
                << "xmid1: " << xmid1 << "\n"
                << "xmid2: " << xmid2 << "\n"
                << "A:\n" << A << "\n"
                << "B:\n" << B << "\n"
                << "A.determinant: " << A.determinant() << "\n"
                << "C:\n" << C << "\n";
      auto k = 6*(i + j*mnx1);
      mcoeffs[k    ] = C(0);
      mcoeffs[k + 1] = C(1);
      mcoeffs[k + 2] = C(2);
      mcoeffs[k + 3] = C(3);
      mcoeffs[k + 4] = C(4);
      mcoeffs[k + 5] = C(5);
    }
  }
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::operator()(const Vector& pos) const {
  const auto i0 = lowerBound(pos);
  const auto x = pos[0];
  const auto y = pos[1];
  return mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x + mcoeffs[i0 + 5]*y*y;
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
BiQuadraticInterpolator::lowerBound(const Vector& pos) const {
  const auto x = mxlog ? log(pos[0]) : pos[0];
  const auto y = mylog ? log(pos[1]) : pos[1];
  const auto result = 6u*(mnx1*std::min(mny1, size_t(std::max(0.0, y - mxmin[1])/mxstep[1])) +
                               std::min(mnx1, size_t(std::max(0.0, x - mxmin[0])/mxstep[0])));
  ENSURE(result <= 6u*mnx1*mny1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
BiQuadraticInterpolator::size() const {
  return mcoeffs.size();
}

inline
typename BiQuadraticInterpolator::Vector
BiQuadraticInterpolator::xmin() const {
  return mxmin;
}

inline
typename BiQuadraticInterpolator::Vector
BiQuadraticInterpolator::xmax() const {
  return mxmax;
}

inline
typename BiQuadraticInterpolator::Vector
BiQuadraticInterpolator::xstep() const {
  return mxstep;
}

inline
bool
BiQuadraticInterpolator::xlog() const {
  return mxlog;
}

inline
bool
BiQuadraticInterpolator::ylog() const {
  return mylog;
}

inline
const std::vector<double>&
BiQuadraticInterpolator::coeffs() const {
  return mcoeffs;
}

}
