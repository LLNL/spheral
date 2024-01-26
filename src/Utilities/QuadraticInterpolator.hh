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

#include "QuadraticInterpolatorView.hh"

namespace Spheral {

class QuadraticInterpolator : public QuadraticInterpolatorView{
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F);
  QuadraticInterpolator() { mcoeffs = CoeffsType(0); };
  ~QuadraticInterpolator() {};

  //QuadraticInterpolator(QuadraticInterpolator const& rhs) = default;
  QuadraticInterpolator(QuadraticInterpolator const& rhs) : QuadraticInterpolatorView(rhs) { mcoeffs = deepCopy(rhs.mcoeffs); }

  QuadraticInterpolator& operator=(QuadraticInterpolator const& rhs) {
    if (this != &rhs) {
      mN1 = rhs.mN1;
      mXmin = rhs.mXmin;
      mXmax = rhs.mXmax;
      mcoeffs = deepCopy(rhs.mcoeffs);
    }
    return *this;
  }
  
  bool operator==(QuadraticInterpolator const& rhs) const {
    return ((mN1 == rhs.mN1) and
            (mXmin == rhs.mXmin) and
            (mXmax == rhs.mXmax) and
            (compare(mcoeffs, rhs.mcoeffs)));
  }

  QuadraticInterpolatorView toView() { return QuadraticInterpolatorView(*this); }

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

};

}

#include "QuadraticInterpolatorInline.hh"

#endif
