#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
BiCubicInterpolator::BiCubicInterpolator():
  XYInterpolator() {
}

//------------------------------------------------------------------------------
// Construct by sampling the given functor
//------------------------------------------------------------------------------
template<typename Func>
BiCubicInterpolator::BiCubicInterpolator(const double xmin,
                                         const double xmax,
                                         const double ymin,
                                         const double ymax,
                                         const size_t nx,
                                         const size_t ny,
                                         const Func& F,
                                         const bool xlog,
                                         const bool ylog):
  XYInterpolator(3u, xmin, xmax, ymin, ymax, nx, ny, xlog, ylog) {

  std::function<Dim<2>::SymTensor(const double, const double)> Fgrad = [&](const double x, const double y) -> Dim<2>::SymTensor {
    size_t ix, iy, i0;
    lowerBound(x, y, ix, iy, i0);
    const auto dx = 0.01*(xcoord(ix + 1u) - xcoord(ix));
    const auto dy = 0.01*(ycoord(iy + 1u) - ycoord(iy));
    const auto dFdx_y1 = (F(x + dx, y + dx) - F(x - dx, y + dx))/(2.0*dx);
    const auto dFdx_y0 = (F(x + dx, y - dx) - F(x - dx, y - dx))/(2.0*dx);
    Dim<2>::SymTensor result;
    result.xx((F(x + dx, y) - F(x - dx, y))/(2.0*dx));
    result.yy((F(x, y + dy) - F(x, y - dy))/(2.0*dy));
    result.xy((dFdx_y1 - dFdx_y0)/(2.0*dy));
    return result;
  };

  initialize(F, Fgrad);
}

//   // Pre-determine the A inverse matix
//   double x1 = 0.0, x2 = 1.0/3.0, x3 = 2.0/3.0, x4 = 1.0;
//   double y1 = 0.0, y2 = 1.0/3.0, y3 = 2.0/3.0, y4 = 1.0;
//   Eigen::MatrixXd A(16, 16);
//   A <<
//     1.0, x1, x1*x1, x1*x1*x1, y1, x1*y1, x1*x1*y1, x1*x1*x1*y1, y1*y1, x1*y1*y1, x1*x1*y1*y1, x1*x1*x1*y1*y1, y1*y1*y1, x1*y1*y1*y1, x1*x1*y1*y1*y1, x1*x1*x1*y1*y1*y1,
//     1.0, x2, x2*x2, x2*x2*x2, y1, x2*y1, x2*x2*y1, x2*x2*x2*y1, y1*y1, x2*y1*y1, x2*x2*y1*y1, x2*x2*x2*y1*y1, y1*y1*y1, x2*y1*y1*y1, x2*x2*y1*y1*y1, x2*x2*x2*y1*y1*y1,
//     1.0, x3, x3*x3, x3*x3*x3, y1, x3*y1, x3*x3*y1, x3*x3*x3*y1, y1*y1, x3*y1*y1, x3*x3*y1*y1, x3*x3*x3*y1*y1, y1*y1*y1, x3*y1*y1*y1, x3*x3*y1*y1*y1, x3*x3*x3*y1*y1*y1,
//     1.0, x4, x4*x4, x4*x4*x4, y1, x4*y1, x4*x4*y1, x4*x4*x4*y1, y1*y1, x4*y1*y1, x4*x4*y1*y1, x4*x4*x4*y1*y1, y1*y1*y1, x4*y1*y1*y1, x4*x4*y1*y1*y1, x4*x4*x4*y1*y1*y1,

//     1.0, x1, x1*x1, x1*x1*x1, y2, x1*y2, x1*x1*y2, x1*x1*x1*y2, y2*y2, x1*y2*y2, x1*x1*y2*y2, x1*x1*x1*y2*y2, y2*y2*y2, x1*y2*y2*y2, x1*x1*y2*y2*y2, x1*x1*x1*y2*y2*y2,
//     1.0, x2, x2*x2, x2*x2*x2, y2, x2*y2, x2*x2*y2, x2*x2*x2*y2, y2*y2, x2*y2*y2, x2*x2*y2*y2, x2*x2*x2*y2*y2, y2*y2*y2, x2*y2*y2*y2, x2*x2*y2*y2*y2, x2*x2*x2*y2*y2*y2,
//     1.0, x3, x3*x3, x3*x3*x3, y2, x3*y2, x3*x3*y2, x3*x3*x3*y2, y2*y2, x3*y2*y2, x3*x3*y2*y2, x3*x3*x3*y2*y2, y2*y2*y2, x3*y2*y2*y2, x3*x3*y2*y2*y2, x3*x3*x3*y2*y2*y2,
//     1.0, x4, x4*x4, x4*x4*x4, y2, x4*y2, x4*x4*y2, x4*x4*x4*y2, y2*y2, x4*y2*y2, x4*x4*y2*y2, x4*x4*x4*y2*y2, y2*y2*y2, x4*y2*y2*y2, x4*x4*y2*y2*y2, x4*x4*x4*y2*y2*y2,

//     1.0, x1, x1*x1, x1*x1*x1, y3, x1*y3, x1*x1*y3, x1*x1*x1*y3, y3*y3, x1*y3*y3, x1*x1*y3*y3, x1*x1*x1*y3*y3, y3*y3*y3, x1*y3*y3*y3, x1*x1*y3*y3*y3, x1*x1*x1*y3*y3*y3,
//     1.0, x2, x2*x2, x2*x2*x2, y3, x2*y3, x2*x2*y3, x2*x2*x2*y3, y3*y3, x2*y3*y3, x2*x2*y3*y3, x2*x2*x2*y3*y3, y3*y3*y3, x2*y3*y3*y3, x2*x2*y3*y3*y3, x2*x2*x2*y3*y3*y3,
//     1.0, x3, x3*x3, x3*x3*x3, y3, x3*y3, x3*x3*y3, x3*x3*x3*y3, y3*y3, x3*y3*y3, x3*x3*y3*y3, x3*x3*x3*y3*y3, y3*y3*y3, x3*y3*y3*y3, x3*x3*y3*y3*y3, x3*x3*x3*y3*y3*y3,
//     1.0, x4, x4*x4, x4*x4*x4, y3, x4*y3, x4*x4*y3, x4*x4*x4*y3, y3*y3, x4*y3*y3, x4*x4*y3*y3, x4*x4*x4*y3*y3, y3*y3*y3, x4*y3*y3*y3, x4*x4*y3*y3*y3, x4*x4*x4*y3*y3*y3,

//     1.0, x1, x1*x1, x1*x1*x1, y4, x1*y4, x1*x1*y4, x1*x1*x1*y4, y4*y4, x1*y4*y4, x1*x1*y4*y4, x1*x1*x1*y4*y4, y4*y4*y4, x1*y4*y4*y4, x1*x1*y4*y4*y4, x1*x1*x1*y4*y4*y4,
//     1.0, x2, x2*x2, x2*x2*x2, y4, x2*y4, x2*x2*y4, x2*x2*x2*y4, y4*y4, x2*y4*y4, x2*x2*y4*y4, x2*x2*x2*y4*y4, y4*y4*y4, x2*y4*y4*y4, x2*x2*y4*y4*y4, x2*x2*x2*y4*y4*y4,
//     1.0, x3, x3*x3, x3*x3*x3, y4, x3*y4, x3*x3*y4, x3*x3*x3*y4, y4*y4, x3*y4*y4, x3*x3*y4*y4, x3*x3*x3*y4*y4, y4*y4*y4, x3*y4*y4*y4, x3*x3*y4*y4*y4, x3*x3*x3*y4*y4*y4,
//     1.0, x4, x4*x4, x4*x4*x4, y4, x4*y4, x4*x4*y4, x4*x4*x4*y4, y4*y4, x4*y4*y4, x4*x4*y4*y4, x4*x4*x4*y4*y4, y4*y4*y4, x4*y4*y4*y4, x4*x4*y4*y4*y4, x4*x4*x4*y4*y4*y4;
//   const auto Ainv = A.inverse();

//   // Fit the coefficients
//   double dx, dy;
//   Eigen::VectorXd b(16), c(16);
//   for (auto i = 0u; i < mnx1; ++i) {
//     for (auto j = 0u; j < mny1; ++j) {
//       x1 = xcoord(i);
//       x4 = xcoord(i + 1u);
//       dx = x4 - x1;
//       x2 = x1 + dx/3.0;
//       x3 = x1 + 2.0*dx/3.0;
//       y1 = ycoord(j);
//       y4 = ycoord(j + 1u);
//       dy = y4 - y1;
//       y2 = y1 + dy/3.0;
//       y3 = y1 + 2.0*dy/3.0;
//       b << F(x1, y1),
//            F(x2, y1),
//            F(x3, y1),
//            F(x4, y1),
//            F(x1, y2),
//            F(x2, y2),
//            F(x3, y2),
//            F(x4, y2),
//            F(x1, y3),
//            F(x2, y3),
//            F(x3, y3),
//            F(x4, y3),
//            F(x1, y4),
//            F(x2, y4),
//            F(x3, y4),
//            F(x4, y4);
//       c = Ainv*b;
//       const auto k = mncoeffs*(i + j*mnx1);
//       mcoeffs[k     ] = c(0);
//       mcoeffs[k +  1] = c(1);
//       mcoeffs[k +  2] = c(2);
//       mcoeffs[k +  3] = c(3);
//       mcoeffs[k +  4] = c(4);
//       mcoeffs[k +  5] = c(5);
//       mcoeffs[k +  6] = c(6);
//       mcoeffs[k +  7] = c(7);
//       mcoeffs[k +  8] = c(8);
//       mcoeffs[k +  9] = c(9);
//       mcoeffs[k + 10] = c(10);
//       mcoeffs[k + 11] = c(11);
//       mcoeffs[k + 12] = c(12);
//       mcoeffs[k + 13] = c(13);
//       mcoeffs[k + 14] = c(14);
//       mcoeffs[k + 15] = c(15);
//     }
//   }
// }

//------------------------------------------------------------------------------
// Construct given a functor and its gradient
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
BiCubicInterpolator::BiCubicInterpolator(const double xmin,
                                         const double xmax,
                                         const double ymin,
                                         const double ymax,
                                         const size_t nx,
                                         const size_t ny,
                                         const Func& F,
                                         const GradFunc& Fgrad,
                                         const bool xlog,
                                         const bool ylog):
  XYInterpolator(3u, xmin, xmax, ymin, ymax, nx, ny, xlog, ylog) {
  initialize(F, Fgrad);
}

//------------------------------------------------------------------------------
// The actual initialization method.
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
void
BiCubicInterpolator::initialize(const Func& F,
                                const GradFunc& Fgrad) {

  // Pre-determine the A inverse matix
  double x0 = 0.0, x1 = 1.0, y0 = 0.0, y1 = 1.0;
  Dim<2>::SymTensor Fgrad00, Fgrad01, Fgrad10, Fgrad11;
  Eigen::MatrixXd A(16, 16), Ainv(16, 16);
  // A <<
  //   1.0, x0, x0*x0, x0*x0*x0, y0, x0*y0, x0*x0*y0, x0*x0*x0*y0, y0*y0, x0*y0*y0, x0*x0*y0*y0, x0*x0*x0*y0*y0, y0*y0*y0, x0*y0*y0*y0, x0*x0*y0*y0*y0, x0*x0*x0*y0*y0*y0,  // f(x0,y0)
  //   1.0, x1, x1*x1, x1*x1*x1, y0, x1*y0, x1*x1*y0, x1*x1*x1*y0, y0*y0, x1*y0*y0, x1*x1*y0*y0, x1*x1*x1*y0*y0, y0*y0*y0, x1*y0*y0*y0, x1*x1*y0*y0*y0, x1*x1*x1*y0*y0*y0,  // f(x1,y0)
  //   1.0, x0, x0*x0, x0*x0*x0, y1, x0*y1, x0*x0*y1, x0*x0*x0*y1, y1*y1, x0*y1*y1, x0*x0*y1*y1, x0*x0*x0*y1*y1, y1*y1*y1, x0*y1*y1*y1, x0*x0*y1*y1*y1, x0*x0*x0*y1*y1*y1,  // f(x0,y1)
  //   1.0, x1, x1*x1, x1*x1*x1, y1, x1*y1, x1*x1*y1, x1*x1*x1*y1, y1*y1, x1*y1*y1, x1*x1*y1*y1, x1*x1*x1*y1*y1, y1*y1*y1, x1*y1*y1*y1, x1*x1*y1*y1*y1, x1*x1*x1*y1*y1*y1,  // f(x1,y1)

  //   0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, y0, 2.0*x0*y0, 3.0*x0*x0*y0, 0.0, y0*y0, 2.0*x0*y0*y0, 3.0*x0*x0*y0*y0, 0.0, y0*y0*y0, 2.0*x0*y0*y0*y0, 3.0*x0*x0*y0*y0*y0,        // \partial_x f(x0,y0)
  //   0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, y0, 2.0*x1*y0, 3.0*x1*x1*y0, 0.0, y0*y0, 2.0*x1*y0*y0, 3.0*x1*x1*y0*y0, 0.0, y0*y0*y0, 2.0*x1*y0*y0*y0, 3.0*x1*x1*y0*y0*y0,        // \partial_x f(x1,y0)    
  //   0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, y1, 2.0*x0*y1, 3.0*x0*x0*y1, 0.0, y1*y1, 2.0*x0*y1*y1, 3.0*x0*x0*y1*y1, 0.0, y1*y1*y1, 2.0*x0*y1*y1*y1, 3.0*x0*x0*y1*y1*y1,        // \partial_x f(x0,y1)
  //   0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, y1, 2.0*x1*y1, 3.0*x1*x1*y1, 0.0, y1*y1, 2.0*x1*y1*y1, 3.0*x1*x1*y1*y1, 0.0, y1*y1*y1, 2.0*x1*y1*y1*y1, 3.0*x1*x1*y1*y1*y1,        // \partial_x f(x1,y1)
    
  //   0.0, 0.0, 0.0, 0.0, 1.0, x0, x0*x0, x0*x0*x0, 2.0*y0, 2.0*x0*y0, 2.0*x0*x0*y0, 2.0*x0*x0*x0*y0, 3.0*y0*y0, 3.0*x0*y0*y0, 3.0*x0*x0*y0*y0, 3.0*x0*x0*x0*y0*y0,        // \partial_y f(x0,y0)
  //   0.0, 0.0, 0.0, 0.0, 1.0, x1, x1*x1, x1*x1*x1, 2.0*y0, 2.0*x1*y0, 2.0*x1*x1*y0, 2.0*x1*x1*x1*y0, 3.0*y0*y0, 3.0*x1*y0*y0, 3.0*x1*x1*y0*y0, 3.0*x1*x1*x1*y0*y0,        // \partial_y f(x1,y0)
  //   0.0, 0.0, 0.0, 0.0, 1.0, x0, x0*x0, x0*x0*x0, 2.0*y1, 2.0*x0*y1, 2.0*x0*x0*y1, 2.0*x0*x0*x0*y1, 3.0*y1*y1, 3.0*x0*y1*y1, 3.0*x0*x0*y1*y1, 3.0*x0*x0*x0*y1*y1,        // \partial_y f(x0,y1)
  //   0.0, 0.0, 0.0, 0.0, 1.0, x1, x1*x1, x1*x1*x1, 2.0*y1, 2.0*x1*y1, 2.0*x1*x1*y1, 2.0*x1*x1*x1*y1, 3.0*y1*y1, 3.0*x1*y1*y1, 3.0*x1*x1*y1*y1, 3.0*x1*x1*x1*y1*y1,        // \partial_y f(x1,y1)

  //   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, 2.0*y0, 4.0*x0*y0, 6.0*x0*x0*y0, 0.0, 3.0*y0*y0, 6.0*x0*y0*y0, 9.0*x0*x0*y0*y0,                                // \partial_xy f(x0,y0)
  //   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, 2.0*y0, 4.0*x1*y0, 6.0*x1*x1*y0, 0.0, 3.0*y0*y0, 6.0*x1*y0*y0, 9.0*x1*x1*y0*y0,                                // \partial_xy f(x1,y0)
  //   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, 2.0*y1, 4.0*x0*y1, 6.0*x0*x0*y1, 0.0, 3.0*y1*y1, 6.0*x0*y1*y1, 9.0*x0*x0*y1*y1,                                // \partial_xy f(x0,y1)
  //   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, 2.0*y1, 4.0*x1*y1, 6.0*x1*x1*y1, 0.0, 3.0*y1*y1, 6.0*x1*y1*y1, 9.0*x1*x1*y1*y1;                                // \partial_xy f(x1,y1)
  // const auto Ainv = A.inverse();

  // Eigen::MatrixXd Ainv(16, 16);
  // Ainv <<  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  //          0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  //         -3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  //          2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  
  //          0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  
  //          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  
  //          0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0,  
  //          0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  
  //         -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  
  //          0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  
  //          9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1,
  //         -6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1, 
  //          2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0, 
  //          0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  
  //         -6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1, 
  //          4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1;

  // Fit the coefficients
  Eigen::VectorXd b(16), c(16);
  double dx, dy;
  for (auto i = 0u; i < mnx1; ++i) {
    x0 = xcoord(i);
    x1 = xcoord(i + 1u);
    dx = x1 - x0;
    for (auto j = 0u; j < mny1; ++j) {
      y0 = ycoord(j);
      y1 = ycoord(j + 1u);
      dy = y1 - y0;
      Fgrad00 = Fgrad(x0, y0);
      Fgrad10 = Fgrad(x1, y0);
      Fgrad01 = Fgrad(x0, y1);
      Fgrad11 = Fgrad(x1, y1);
      b << F(x0, y0),
           F(x1, y0),
           F(x0, y1),
           F(x1, y1),
           Fgrad00.xx(),// * dx,     // partial_x
           Fgrad10.xx(),// * dx,     // partial_x
           Fgrad01.xx(),// * dx,     // partial_x
           Fgrad11.xx(),// * dx,     // partial_x
           Fgrad00.yy(),// * dy,     // partial_y
           Fgrad10.yy(),// * dy,     // partial_y
           Fgrad01.yy(),// * dy,     // partial_y
           Fgrad11.yy(),// * dy,     // partial_y
           Fgrad00.xy(),// * dx*dy,  // partial_xy
           Fgrad10.xy(),// * dx*dy,  // partial_xy
           Fgrad01.xy(),// * dx*dy,  // partial_xy
           Fgrad11.xy();// * dx*dy;  // partial_xy
      CHECK(b == b);  // Looking for NaNs

      A <<
        1.0, x0, x0*x0, x0*x0*x0, y0, x0*y0, x0*x0*y0, x0*x0*x0*y0, y0*y0, x0*y0*y0, x0*x0*y0*y0, x0*x0*x0*y0*y0, y0*y0*y0, x0*y0*y0*y0, x0*x0*y0*y0*y0, x0*x0*x0*y0*y0*y0,  // f(x0,y0)
        1.0, x1, x1*x1, x1*x1*x1, y0, x1*y0, x1*x1*y0, x1*x1*x1*y0, y0*y0, x1*y0*y0, x1*x1*y0*y0, x1*x1*x1*y0*y0, y0*y0*y0, x1*y0*y0*y0, x1*x1*y0*y0*y0, x1*x1*x1*y0*y0*y0,  // f(x1,y0)
        1.0, x0, x0*x0, x0*x0*x0, y1, x0*y1, x0*x0*y1, x0*x0*x0*y1, y1*y1, x0*y1*y1, x0*x0*y1*y1, x0*x0*x0*y1*y1, y1*y1*y1, x0*y1*y1*y1, x0*x0*y1*y1*y1, x0*x0*x0*y1*y1*y1,  // f(x0,y1)
        1.0, x1, x1*x1, x1*x1*x1, y1, x1*y1, x1*x1*y1, x1*x1*x1*y1, y1*y1, x1*y1*y1, x1*x1*y1*y1, x1*x1*x1*y1*y1, y1*y1*y1, x1*y1*y1*y1, x1*x1*y1*y1*y1, x1*x1*x1*y1*y1*y1,  // f(x1,y1)

        0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, y0, 2.0*x0*y0, 3.0*x0*x0*y0, 0.0, y0*y0, 2.0*x0*y0*y0, 3.0*x0*x0*y0*y0, 0.0, y0*y0*y0, 2.0*x0*y0*y0*y0, 3.0*x0*x0*y0*y0*y0,        // \partial_x f(x0,y0)
        0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, y0, 2.0*x1*y0, 3.0*x1*x1*y0, 0.0, y0*y0, 2.0*x1*y0*y0, 3.0*x1*x1*y0*y0, 0.0, y0*y0*y0, 2.0*x1*y0*y0*y0, 3.0*x1*x1*y0*y0*y0,        // \partial_x f(x1,y0)    
        0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, y1, 2.0*x0*y1, 3.0*x0*x0*y1, 0.0, y1*y1, 2.0*x0*y1*y1, 3.0*x0*x0*y1*y1, 0.0, y1*y1*y1, 2.0*x0*y1*y1*y1, 3.0*x0*x0*y1*y1*y1,        // \partial_x f(x0,y1)
        0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, y1, 2.0*x1*y1, 3.0*x1*x1*y1, 0.0, y1*y1, 2.0*x1*y1*y1, 3.0*x1*x1*y1*y1, 0.0, y1*y1*y1, 2.0*x1*y1*y1*y1, 3.0*x1*x1*y1*y1*y1,        // \partial_x f(x1,y1)
    
        0.0, 0.0, 0.0, 0.0, 1.0, x0, x0*x0, x0*x0*x0, 2.0*y0, 2.0*x0*y0, 2.0*x0*x0*y0, 2.0*x0*x0*x0*y0, 3.0*y0*y0, 3.0*x0*y0*y0, 3.0*x0*x0*y0*y0, 3.0*x0*x0*x0*y0*y0,        // \partial_y f(x0,y0)
        0.0, 0.0, 0.0, 0.0, 1.0, x1, x1*x1, x1*x1*x1, 2.0*y0, 2.0*x1*y0, 2.0*x1*x1*y0, 2.0*x1*x1*x1*y0, 3.0*y0*y0, 3.0*x1*y0*y0, 3.0*x1*x1*y0*y0, 3.0*x1*x1*x1*y0*y0,        // \partial_y f(x1,y0)
        0.0, 0.0, 0.0, 0.0, 1.0, x0, x0*x0, x0*x0*x0, 2.0*y1, 2.0*x0*y1, 2.0*x0*x0*y1, 2.0*x0*x0*x0*y1, 3.0*y1*y1, 3.0*x0*y1*y1, 3.0*x0*x0*y1*y1, 3.0*x0*x0*x0*y1*y1,        // \partial_y f(x0,y1)
        0.0, 0.0, 0.0, 0.0, 1.0, x1, x1*x1, x1*x1*x1, 2.0*y1, 2.0*x1*y1, 2.0*x1*x1*y1, 2.0*x1*x1*x1*y1, 3.0*y1*y1, 3.0*x1*y1*y1, 3.0*x1*x1*y1*y1, 3.0*x1*x1*x1*y1*y1,        // \partial_y f(x1,y1)

        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, 2.0*y0, 4.0*x0*y0, 6.0*x0*x0*y0, 0.0, 3.0*y0*y0, 6.0*x0*y0*y0, 9.0*x0*x0*y0*y0,                                // \partial_xy f(x0,y0)
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, 2.0*y0, 4.0*x1*y0, 6.0*x1*x1*y0, 0.0, 3.0*y0*y0, 6.0*x1*y0*y0, 9.0*x1*x1*y0*y0,                                // \partial_xy f(x1,y0)
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x0, 3.0*x0*x0, 0.0, 2.0*y1, 4.0*x0*y1, 6.0*x0*x0*y1, 0.0, 3.0*y1*y1, 6.0*x0*y1*y1, 9.0*x0*x0*y1*y1,                                // \partial_xy f(x0,y1)
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0*x1, 3.0*x1*x1, 0.0, 2.0*y1, 4.0*x1*y1, 6.0*x1*x1*y1, 0.0, 3.0*y1*y1, 6.0*x1*y1*y1, 9.0*x1*x1*y1*y1;                                // \partial_xy f(x1,y1)


      // c = A.colPivHouseholderQr().solve(b);
      // c = Ainv*b;
      Ainv = A.inverse();
      c = Ainv*b;
      // std::cerr << "BiCubicInterpolator: "
      //           << "\n  F(" << x0 << " " << y0 << ")=" << F(x0, y0)
      //           << "\n  F(" << x1 << " " << y0 << ")=" << F(x1, y0)
      //           << "\n  F(" << x0 << " " << y1 << ")=" << F(x0, y1)
      //           << "\n  F(" << x1 << " " << y1 << ")=" << F(x1, y1)
      //           << "\n  b=" << b << "\n  c=" << c << "\n  A*c=" << A*c << std::endl;
      const auto k = mncoeffs*(i + j*mnx1);
      mcoeffs[k     ] = c(0);
      mcoeffs[k +  1] = c(1);
      mcoeffs[k +  2] = c(2);
      mcoeffs[k +  3] = c(3);
      mcoeffs[k +  4] = c(4);
      mcoeffs[k +  5] = c(5);
      mcoeffs[k +  6] = c(6);
      mcoeffs[k +  7] = c(7);
      mcoeffs[k +  8] = c(8);
      mcoeffs[k +  9] = c(9);
      mcoeffs[k + 10] = c(10);
      mcoeffs[k + 11] = c(11);
      mcoeffs[k + 12] = c(12);
      mcoeffs[k + 13] = c(13);
      mcoeffs[k + 14] = c(14);
      mcoeffs[k + 15] = c(15);
    }
  }
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
BiCubicInterpolator::~BiCubicInterpolator() {
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::operator()(const double xi,
                                const double yi) const {
  double x = xi, y = yi;
  size_t ix, iy, i0;
  lowerBound(xi, yi, ix, iy, i0);
  // eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x, x3 = x*x*x;
  const auto y2 = y*y, y3 = y*y*y;
  return (mcoeffs[i0]               +          // c00
          mcoeffs[i0 +  1]*x        +          // c10
          mcoeffs[i0 +  2]*x2       +          // c20
          mcoeffs[i0 +  3]*x3       +          // c30
          mcoeffs[i0 +  4]    * y   +          // c01
          mcoeffs[i0 +  5]*x  * y   +          // c11
          mcoeffs[i0 +  6]*x2 * y   +          // c21
          mcoeffs[i0 +  7]*x3 * y   +          // c31
          mcoeffs[i0 +  8]    * y2  +          // c02
          mcoeffs[i0 +  9]*x  * y2  +          // c12
          mcoeffs[i0 + 10]*x2 * y2  +          // c22
          mcoeffs[i0 + 11]*x3 * y2  +          // c32
          mcoeffs[i0 + 12]    * y3  +          // c03
          mcoeffs[i0 + 13]*x  * y3  +          // c13
          mcoeffs[i0 + 14]*x2 * y3  +          // c23
          mcoeffs[i0 + 15]*x3 * y3);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime_x(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x;
  const auto y2 = y*y, y3 = y*y*y;
  return (mcoeffs[i0 +  1]              +          // c10
          mcoeffs[i0 +  2]*2.0*x        +          // c20
          mcoeffs[i0 +  3]*3.0*x2       +          // c30
          mcoeffs[i0 +  5]        * y   +          // c11
          mcoeffs[i0 +  6]*2.0*x  * y   +          // c21
          mcoeffs[i0 +  7]*3.0*x2 * y   +          // c31
          mcoeffs[i0 +  9]        * y2  +          // c12
          mcoeffs[i0 + 10]*2.0*x  * y2  +          // c22
          mcoeffs[i0 + 11]*3.0*x2 * y2  +          // c32
          mcoeffs[i0 + 13]        * y3  +          // c13
          mcoeffs[i0 + 14]*2.0*x  * y3  +          // c23
          mcoeffs[i0 + 15]*3.0*x2 * y3);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime_y(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x, x3 = x*x*x;
  const auto y2 = y*y;
  return (mcoeffs[i0 +  4]              +          // c01
          mcoeffs[i0 +  5]*x            +          // c11
          mcoeffs[i0 +  6]*x2           +          // c21
          mcoeffs[i0 +  7]*x3           +          // c31
          mcoeffs[i0 +  8]    * 2.0*y   +          // c02
          mcoeffs[i0 +  9]*x  * 2.0*y   +          // c12
          mcoeffs[i0 + 10]*x2 * 2.0*y   +          // c22
          mcoeffs[i0 + 11]*x3 * 2.0*y   +          // c32
          mcoeffs[i0 + 12]    * 3.0*y2  +          // c03
          mcoeffs[i0 + 13]*x  * 3.0*y2  +          // c13
          mcoeffs[i0 + 14]*x2 * 3.0*y2  +          // c23
          mcoeffs[i0 + 15]*x3 * 3.0*y2);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xx)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_xx(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto y2 = y*y, y3 = y*y*y;
  return (mcoeffs[i0 +  2]*2.0         +          // c20
          mcoeffs[i0 +  3]*6.0*x       +          // c30
          mcoeffs[i0 +  6]*2.0         +          // c21
          mcoeffs[i0 +  7]*6.0*x * y   +          // c31
          mcoeffs[i0 + 10]*2.0   * y2  +          // c22
          mcoeffs[i0 + 11]*6.0*x * y2  +          // c32
          mcoeffs[i0 + 14]*2.0   * y3  +          // c23
          mcoeffs[i0 + 15]*6.0*x * y3);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xy)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_xy(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x;
  const auto y2 = y*y;
  return (mcoeffs[i0 +  5]   +                         // c11
          mcoeffs[i0 +  6]*2.0*x            +          // c21
          mcoeffs[i0 +  7]*3.0*x            +          // c31
          mcoeffs[i0 +  9]        * 2.0*y   +          // c12
          mcoeffs[i0 + 10]*2.0*x  * 2.0*y   +          // c22
          mcoeffs[i0 + 11]*3.0*x2 * 2.0*y   +          // c32
          mcoeffs[i0 + 13]        * 3.0*y2  +          // c13
          mcoeffs[i0 + 14]*2.0*x  * 3.0*y2  +          // c23
          mcoeffs[i0 + 15]*3.0*x2 * 3.0*y2);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (yy)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_yy(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x, x3 = x*x*x;
  return (mcoeffs[i0 +  8]    * 2.0    +          // c02
          mcoeffs[i0 +  9]*x  * 2.0    +          // c12
          mcoeffs[i0 + 10]*x2 * 2.0    +          // c22
          mcoeffs[i0 + 11]*x3 * 2.0    +          // c32
          mcoeffs[i0 + 12]    * 6.0*y  +          // c03
          mcoeffs[i0 + 13]*x  * 6.0*y  +          // c13
          mcoeffs[i0 + 14]*x2 * 6.0*y  +          // c23
          mcoeffs[i0 + 15]*x3 * 6.0*y);           // c33
}

}
