//---------------------------------Spheral++----------------------------------//
// NSincPolynomialKernel -- The sinc interpolation kernel: W = sin(pi*eta)/(pi*eta).
//
// Created by JMO, Tue Jan  7 15:01:13 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_NSincPolynomialKernel_hh__
#define __Spheral_NSincPolynomialKernel_hh__

#include "Kernel.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class NSincPolynomialKernel: 
    public Kernel<Dimension, NSincPolynomialKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  NSincPolynomialKernel(const int order);
  ~NSincPolynomialKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaMagnitude, const double Hdet) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Order of the polynomials.
  int mOrder;

#ifndef __GCCXML__
  // The coefficients for each piecewise section.
  std::vector< std::vector<double> > mAij;

  // Private method to fill in the polynomial coefficients.
  void setPolynomialCoefficients(const int order, 
                                 std::vector< std::vector<double> >& Aij) const;
#endif

};

}

#include "NSincPolynomialKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class NSincPolynomialKernel;
}

#endif
