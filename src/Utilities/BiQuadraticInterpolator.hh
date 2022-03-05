//---------------------------------Spheral++----------------------------------//
// BiQuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 2D
// Assumes the results is interpolated as
//   <F(x,y)> = c0 + c1*x + c2*y + c3*x*y + c4*x^2 + c5*x^2*y + c6*y^2 + c7*x*y^2 + c8*x^2*y^2
//
// Created by JMO, Thu Dec 10 14:48:01 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_BiQuadraticInterpolator__
#define __Spheral_BiQuadraticInterpolator__

#include <vector>

namespace Spheral {

class BiQuadraticInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  BiQuadraticInterpolator(const double xmin,
                          const double xmax,
                          const double ymin,
                          const double ymax,
                          const size_t nx,
                          const size_t ny,
                          const Func& F,
                          const bool xlog = false,
                          const bool ylog = false);
  BiQuadraticInterpolator();
  ~BiQuadraticInterpolator();

  // Initialize for interpolating the given function
  template<typename Func>
  void initialize(const double xmin,
                  const double xmax,
                  const double ymin,
                  const double ymax,
                  const size_t nx,
                  const size_t ny,
                  const Func& F,
                  const bool xlog = false,
                  const bool ylog = false);

  // Interpolate for the F(x,y) value
  double operator()(const double x, const double y) const;

  // Interpolated gradient values.
  double prime_x(const double x, const double y) const;
  double prime_y(const double x, const double y) const;
  double prime2_xx(const double x, const double y) const;
  double prime2_xy(const double x, const double y) const;
  double prime2_yx(const double x, const double y) const;
  double prime2_yy(const double x, const double y) const;

  // Return the lower bound index in the table of coefficients for the given position
  size_t lowerBound(const double x, const double y) const;

  // Allow read access the internal data representation
  size_t size() const;                        // The size of the tabulated coefficient arrays
  double xmin() const;                        // Minimum x coordinate for table              
  double xmax() const;                        // Maximum x coordinate for table              
  double ymin() const;                        // Minimum y coordinate for table              
  double ymax() const;                        // Maximum y coordinate for table              
  double xstep() const;                       // x step size
  double ystep() const;                       // y step size
  bool xlog() const;                          // Are we using log spacing in x?
  bool ylog() const;                          // Are we using log spacing in y?
  const std::vector<double>& coeffs() const;  // the fitting coefficients
  
  // Comparison
  bool operator==(const BiQuadraticInterpolator& rhs) const;

private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  bool mxlog, mylog;
  size_t mnx1, mny1;
  double mxmin, mxmax, mymin, mymax, mxstep, mystep;
  std::vector<double> mcoeffs;

  // Compute a coordinate value depending on whether we're using log-space
  double coord(const double xmin, const double dx,
               const size_t ix, const size_t nx,
               const bool xlog) const;
};

}

#include "BiQuadraticInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class BiQuadraticInterpolator;
}

#endif
