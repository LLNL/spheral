//---------------------------------Spheral++----------------------------------//
// XYInterpolator
//
// Provides utility functions for use with our 2D spatial interpolators
//
// Created by JMO, Sun Mar 13 21:31:09 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_XYInterpolator__
#define __Spheral_XYInterpolator__

#include <vector>

namespace Spheral {

class XYInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  XYInterpolator(const size_t order,
                 const double xmin,
                 const double xmax,
                 const double ymin,
                 const double ymax,
                 const size_t nx,
                 const size_t ny,
                 const bool xlog,
                 const bool ylog);
                 
  XYInterpolator();
  virtual ~XYInterpolator();

  // Return the lower bound index in the table of coefficients for the given position
  void lowerBound(const double x, const double y,
                  size_t& ix, size_t& iy, size_t& i0) const;

  // Allow read access the internal data representation
  double xmin() const;                        // Minimum x coordinate for table
  double xmax() const;                        // Maximum x coordinate for table
  double ymin() const;                        // Minimum y coordinate for table
  double ymax() const;                        // Maximum y coordinate for table
  double xstep() const;                       // x step size
  double ystep() const;                       // y step size
  bool xlog() const;                          // Are we using log spacing in x?
  bool ylog() const;                          // Are we using log spacing in y?
  size_t size() const;                        // The size of the tabulated coefficient arrays
  const std::vector<double>& coeffs() const;  // the fitting coefficients
  
  // Comparison
  bool operator==(const XYInterpolator& rhs) const;

protected:
  //--------------------------- Protected Interface --------------------------//
  // Member data
  bool mxlog, mylog;
  size_t mnx1, mny1, mncoeffs;
  double mxmin, mxmax, mymin, mymax, mxstep, mystep;
  std::vector<double> mcoeffs;

  // Compute a coordinate value depending on whether we're using log-space
  double coord(const double xmin, const double dx,
               const size_t ix, const size_t nx,
               const bool xlog) const;
  double xcoord(const size_t ix) const;
  double ycoord(const size_t iy) const;

  // Similar to above, but compute the relative normalized coordinate inside
  // a grid patch for the fit (range [0,1]).
  void eta_coords(const double xi, const double yi,
                  double& etax,
                  double& etay,
                  size_t& ix,
                  size_t& iy,
                  size_t& i0) const;
};

}

#include "XYInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class XYInterpolator;
}

#endif
