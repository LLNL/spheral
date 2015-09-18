//---------------------------------Spheral++----------------------------------//
// QuinticSplineKernel -- A quintic spline, as described in
// Bonet & Kulasegaruam 2002, Appl. Math. Comput., 126, 133-155.
//
// Kernel extent: 2.0
//
// Created by JMO, Wed Jul  9 16:24:25 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuinticSplineKernel_hh__
#define __Spheral_QuinticSplineKernel_hh__

#include "Kernel.hh"

namespace Spheral {
namespace KernelSpace {

template<typename Dimension>
class QuinticSplineKernel: public Kernel<Dimension, QuinticSplineKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  QuinticSplineKernel();
  ~QuinticSplineKernel();

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
#include "QuinticSplineKernelInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class QuinticSplineKernel;
  }
}

#endif
