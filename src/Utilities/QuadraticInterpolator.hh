//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuadraticInterpolator__
#define __Spheral_QuadraticInterpolator__

#include "chai/ManagedArray.hpp"
#include "config.hh"

#include <cstddef>
#include <vector>

namespace Spheral {

class QIBase {
public:
  using ContainerType = typename chai::ManagedArray<double>;
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  SPHERAL_HOST_DEVICE QIBase() = default;

  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const QIBase& rhs) const;

  // Interpolate for the y value
  SPHERAL_HOST_DEVICE double operator()(const double x) const;
  SPHERAL_HOST_DEVICE double prime(const double x) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x) const;   // Second derivative
  // Index access
  SPHERAL_HOST_DEVICE double operator[](const size_t i) const;

  // Same as above, but use a pre-computed table position (from lowerBound)
  SPHERAL_HOST_DEVICE double operator()(const double x, const size_t i0) const;
  SPHERAL_HOST_DEVICE double prime(const double x, const size_t i0) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x, const size_t i0) const;   // Second derivative

  // Return the lower bound index in the table for the given x coordinate
  SPHERAL_HOST_DEVICE size_t lowerBound(const double x) const;

  // Allow read access the internal data representation
  SPHERAL_HOST_DEVICE size_t size() const;                        // The size of the tabulated coefficient arrays
  SPHERAL_HOST_DEVICE double xmin() const;                        // Minimum x coordinate for table
  SPHERAL_HOST_DEVICE double xmax() const;                        // Maximum x coordinate for table
  SPHERAL_HOST_DEVICE double xstep() const;                       // delta x between tabulated values
  void move(chai::ExecutionSpace space) { mcoeffs.move(space); }
  SPHERAL_HOST_DEVICE double* data() const { return mcoeffs.data(); }

  SPHERAL_HOST QIBase(size_t N1,
                      double xmin,
                      double xmax,
                      double xstep,
                      ContainerType const& vals) :
    mN1(N1),
    mXmin(xmin),
    mXmax(xmax),
    mXstep(xstep),
    mcoeffs(vals) { mcoeffs.registerTouch(chai::CPU); }
protected:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1 = 0u;
  double mXmin = 0.;
  double mXmax = 0.;
  double mXstep = 0.;
  ContainerType mcoeffs;
};

class QuadraticInterpolator : public QIBase {
public:
  using ContainerType = typename chai::ManagedArray<double>;
  template<typename Func>
  QuadraticInterpolator(double xmin, double xmax, size_t n, const Func& F);
  QuadraticInterpolator(double xmin, double xmax, const std::vector<double>& yvals);
  QuadraticInterpolator() = default;
  ~QuadraticInterpolator();
  QuadraticInterpolator(const QuadraticInterpolator& rhs);
  QuadraticInterpolator& operator=(const QuadraticInterpolator& rhs);

  // Initialize after construction, either with a function or tabulated values
  template<typename Func>
  void initialize(double xmin, double xmax, size_t n, const Func& f);
  void initialize(double xmin, double xmax, const std::vector<double>& yvals);

  QIBase view() {
    return QIBase(mN1, mXmin, mXmax, mXstep, mcoeffs);
  }
private:
  std::vector<double> mVec;
  void initializeMA();
};
}

#include "QuadraticInterpolatorInline.hh"

#endif
