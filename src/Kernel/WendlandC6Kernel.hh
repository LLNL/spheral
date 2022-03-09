//---------------------------------Spheral++----------------------------------//
// WendlandC6Kernel -- .
//
//
// Volume normalizations:
// 1-D:  A = 99/74
// 2-D:  A = 6435/(3856*Pi)
// 3-D:  A = 45045/(31456*Pi)
//
// Kernel extent: 2.0
//
// Created by CDR, Nov 19 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_WendlandC6Kernel_hh__
#define __Spheral_WendlandC6Kernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class WendlandC6Kernel: public Kernel<Dimension, WendlandC6Kernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  WendlandC6Kernel();
  ~WendlandC6Kernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, const double Hdet) const;

  // Return the second derivative value for a given normalized distance or
  // position.
  double grad2Value(double etaMagnitude, const double Hdet) const;

};

}

#include "WendlandC6KernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class WendlandC6Kernel;
}

#endif
