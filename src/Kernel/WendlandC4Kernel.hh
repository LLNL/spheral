//---------------------------------Spheral++----------------------------------//
// WendlandC4Kernel -- .
//
//
// Volume normalizations:
// 1-D:  A = 54/37
// 2-D:  A = 9/(5*Pi)
// 3-D:  A = 165/(112*Pi)
//
// Kernel extent: 2.0
//
// Created by CDR, Nov 5 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_BSplineKernel_hh__
#define __Spheral_BSplineKernel_hh__

#include "Kernel.hh"

namespace Spheral {
namespace KernelSpace {

template<typename Dimension>
class WendlandC4Kernel: public Kernel<Dimension, WendlandC4Kernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  WendlandC4Kernel();
  ~WendlandC4Kernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, double Hdet) const;

  // Return the second derivative value for a given normalized distance or
  // position.
  double grad2Value(double etaMagnitude, double Hdet) const;

};

}
}

#ifndef __GCCXML__
#include "WendlandC4KernelInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class WendlandC4Kernel;
  }
}

#endif
