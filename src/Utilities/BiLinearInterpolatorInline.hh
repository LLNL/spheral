#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
BiLinearInterpolator::BiLinearInterpolator():
  XYInterpolator(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Construct by sampling the given functor
//------------------------------------------------------------------------------
template<typename Func>
BiLinearInterpolator::BiLinearInterpolator(const double xmin,
                                           const double xmax,
                                           const double ymin,
                                           const double ymax,
                                           const size_t nx,
                                           const size_t ny,
                                           const Func& F,
                                           const bool xlog,
                                           const bool ylog):
  XYInterpolator(xmin, xmax, ymin, ymax, nx, ny, xlog, ylog),
  mcoeffs() {

  // Size stuff up.
  mcoeffs.resize(4*mnx1*mny1);

  // We can predetermine A based on a unit square
  double x1 = 0.0, x2 = 1.0;
  double y1 = 0.0, y2 = 1.0;
  Eigen::MatrixXd A(4, 4);
  A << 1.0, x1, y1, x1*y1,
       1.0, x1, y2, x1*y2, 
       1.0, x2, y1, x2*y1,
       1.0, x2, y2, x2*y2;
  const auto Ainv = A.inverse();

  // Fit the coefficients
  Eigen::VectorXd b(4), c(4);
  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x1 = coord(xmin, mxstep, i,      mnx1, xlog);
      x2 = coord(xmin, mxstep, i + 1u, mnx1, xlog);
      y1 = coord(ymin, mystep, j,      mny1, ylog);
      y2 = coord(ymin, mystep, j + 1u, mny1, ylog);
      c << F(x1, y1),
           F(x1, y2),
           F(x2, y1),
           F(x2, y2);
      b = Ainv * c;
      const auto k = 4u*(i + j*mnx1);
      mcoeffs[k    ] = b(0);
      mcoeffs[k + 1] = b(1);
      mcoeffs[k + 2] = b(2);
      mcoeffs[k + 3] = b(3);
    }
  }
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
BiLinearInterpolator::~BiLinearInterpolator() {
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::operator()(const double xi,
                                 const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  // std::cerr << "================================================================================\n"
  //           << "mxlog, mylog : " << mxlog << " " << mylog << "\n"
  //           << "pos   : " << pos << "\n"
  //           << "(x,y) : " << x << " " << y << "\n"
  //           << "i0    : " << i0 << "\n"
  //           << "coeffs: " << mcoeffs[i0] << " " << mcoeffs[i0 + 1] << " " << mcoeffs[i0 + 2] << " " << mcoeffs[i0 + 3] << " " << mcoeffs[i0 + 4] << " " << mcoeffs[i0 + 5] << "\n"
  //           << "F(x,y): " << mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x + mcoeffs[i0 + 5]*y*y << "\n";
  return (mcoeffs[i0] +
          mcoeffs[i0 + 1]*x +
          mcoeffs[i0 + 2]*y +
          mcoeffs[i0 + 3]*x*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::prime_x(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return (mcoeffs[i0 + 1] +
          mcoeffs[i0 + 3]*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiLinearInterpolator::prime_y(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return (mcoeffs[i0 + 2] +
          mcoeffs[i0 + 3]*x);
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
  return (XYInterpolator::operator==(rhs) and
          (mcoeffs == rhs.mcoeffs));
}

}
