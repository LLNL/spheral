#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
BiLinearInterpolator::BiLinearInterpolator():
  mnx1(),
  mny1(),
  mxmin(),
  mxmax(),
  mymin(),
  mymax(),
  mxstep(),
  mystep(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Construct with tabulated data
//------------------------------------------------------------------------------
template<typename Func>
BiLinearInterpolator::BiLinearInterpolator(const double xmin,
                                                 const double xmax,
                                                 const double ymin,
                                                 const double ymax,
                                                 const size_t nx,
                                                 const size_t ny,
                                                 const Func& F):
  mnx1(),
  mny1(),
  mxmin(),
  mxmax(),
  mxstep(),
  mcoeffs() {
  this->initialize(xmin, xmax, ymin, ymax, nx, ny, F);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
BiLinearInterpolator::~BiLinearInterpolator() {
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
template<typename Func>
inline
void
BiLinearInterpolator::initialize(const double xmin,
                                 const double xmax,
                                 const double ymin,
                                 const double ymax,
                                 const size_t nx,
                                 const size_t ny,
                                 const Func& F) {

  // Size stuff up.
  REQUIRE(nx > 1u);
  REQUIRE(ny > 1u);
  mnx1 = nx - 1u;
  mny1 = ny - 1u;
  mcoeffs.resize(4*mnx1*mny1);

  // Figure out the sampling steps.
  mxmin = xmin;
  mxmax = xmax;
  mymin = ymin;
  mymax = ymax;
  mxstep = (xmax - xmin)/(nx - 1u);
  mystep = (ymax - ymin)/(ny - 1u);

  // Fit the coefficients
  double x1, x2, y1, y2;
  Eigen::MatrixXd A(4, 4);
  Eigen::VectorXd b(4), c(4);
  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x1 = xmin + i      *mxstep;
      x2 = xmin + (i + 1)*mxstep;
      y1 = ymin + j      *mystep;
      y2 = ymin + (j + 1)*mystep;
      A << 1.0, x1, y1, x1*y1,
           1.0, x1, y2, x1*y2, 
           1.0, x2, y1, x2*y1,
           1.0, x2, y2, x2*y2;
      c << F(x1, y1),
           F(x1, y2),
           F(x2, y1),
           F(x2, y2);
      b = A.inverse() * c;
      const auto k = 4*(i + j*mnx1);
      mcoeffs[k    ] = b(0);
      mcoeffs[k + 1] = b(1);
      mcoeffs[k + 2] = b(2);
      mcoeffs[k + 3] = b(3);
    }
  }
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::operator()(const double xi,
                                 const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  // std::cerr << "================================================================================\n"
  //           << "mxlog, mylog : " << mxlog << " " << mylog << "\n"
  //           << "pos   : " << pos << "\n"
  //           << "(x,y) : " << x << " " << y << "\n"
  //           << "i0    : " << i0 << "\n"
  //           << "coeffs: " << mcoeffs[i0] << " " << mcoeffs[i0 + 1] << " " << mcoeffs[i0 + 2] << " " << mcoeffs[i0 + 3] << " " << mcoeffs[i0 + 4] << " " << mcoeffs[i0 + 5] << "\n"
  //           << "F(x,y): " << mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x + mcoeffs[i0 + 5]*y*y << "\n";
  return mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::prime_x(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 1] + mcoeffs[i0 + 3]*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::prime_y(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 2] + mcoeffs[i0 + 3]*x;
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
BiLinearInterpolator::lowerBound(const double x, const double y) const {
  const auto result = 4u*(mnx1*std::min(mny1 - 1u, size_t(std::max(0.0, y - mymin)/mystep)) +
                               std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep)));
  ENSURE(result <= 4u*mnx1*mny1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
BiLinearInterpolator::size() const {
  return mcoeffs.size();
}

inline
double
BiLinearInterpolator::xmin() const {
  return mxmin;
}

inline
double
BiLinearInterpolator::xmax() const {
  return mxmax;
}

inline
double
BiLinearInterpolator::ymin() const {
  return mymin;
}

inline
double
BiLinearInterpolator::ymax() const {
  return mymax;
}

inline
double
BiLinearInterpolator::xstep() const {
  return mxstep;
}

inline
double
BiLinearInterpolator::ystep() const {
  return mystep;
}

inline
const std::vector<double>&
BiLinearInterpolator::coeffs() const {
  return mcoeffs;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
BiLinearInterpolator::operator==(const BiLinearInterpolator& rhs) const {
  return ((mnx1 == rhs.mnx1) and
          (mny1 == rhs.mny1) and
          (mxmin == rhs.mxmin) and
          (mxmax == rhs.mxmax) and
          (mymin == rhs.mymin) and
          (mymax == rhs.mymax) and
          (mxstep == rhs.mxstep) and
          (mystep == rhs.mystep) and
          (mcoeffs == rhs.mcoeffs));
}

}
