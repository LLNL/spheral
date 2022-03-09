//---------------------------------Spheral++----------------------------------//
// ExpInvKernel -- A specialized kernel of the form W = A exp(1/(x + b)).
//                 We currently hardwire b=0.1 and kernelExtent of 3.
//
// Volume normalizations:
// 1-D:  A = 588.13226316830821361 / h
// 2-D:  A = 28.256642009101042845 pi / h^2
// 3-D:  A = 20.151784602855109085 pi / h^3
//
// Kernel extent: 3.0
//
// Created by JMO, Fri Sep 18 14:24:00 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_ExpInvKernel_hh__
#define __Spheral_ExpInvKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class ExpInvKernel: public Kernel<Dimension, ExpInvKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  ExpInvKernel();
  ~ExpInvKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaMagnitude, const double Hdet) const;
};

}

#include "ExpInvKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ExpInvKernel;
}

#endif
