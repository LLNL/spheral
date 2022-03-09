//---------------------------------Spheral++----------------------------------//
// BSplineKernel -- The B spline interpolation kernel.
//
// Monaghan 1992, ARAA, 30, 543
// Monaghan & Lttanzio 1985, A&A, 149, 135
//
// Volume normalizations:
// 1-D:  A = 2/3
// 2-D:  A = 10/(7*Pi)
// 3-D:  A = 1/Pi
//
// Kernel extent: 2.0
//
// Created by JMO, Mon Nov 29 21:51:58 PST 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_BSplineKernel_hh__
#define __Spheral_BSplineKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class BSplineKernel: public Kernel<Dimension, BSplineKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  BSplineKernel();
  ~BSplineKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, const double Hdet) const;

  // Return the second derivative value for a given normalized distance or
  // position.
  double grad2Value(double etaMagnitude, const double Hdet) const;

};

}

#include "BSplineKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class BSplineKernel;
}

#endif
