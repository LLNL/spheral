//---------------------------------Spheral++----------------------------------//
// SphericalKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalKernel_hh__
#define __Spheral_SphericalKernel_hh__

#include "TableKernel.hh"
#include "Utilities/BiQuadraticInterpolator.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalKernel {

public:
  //--------------------------- Public Interface ---------------------------//
  using InterpolatorType = BiQuadraticInterpolator;
  using Scalar = Dim<1>::Scalar;
  using Vector = Dim<1>::Vector;
  using Tensor = Dim<1>::Tensor;
  using SymTensor = Dim<1>::SymTensor;

  // Constructor.
  // Takes a normal 3D TableKernel and constructs the integral form appropriate
  // for 1D spherical coordinates.
  template<typename KernelType> SphericalKernel(const KernelType& kernel,
                                                const unsigned numIntegral = 1000u,
                                                const unsigned numKernel = 200u);
  SphericalKernel(const SphericalKernel& rhs);

  // Take an arbitrary Kernel

  // Destructor.
  virtual ~SphericalKernel();

  // Assignment.
  SphericalKernel& operator=(const SphericalKernel& rhs);

  // Comparisons
  bool operator==(const SphericalKernel& rhs) const;

  // These methods taking a Vector eta and Vector position are the special methods
  // allowing this kernel to implement the asymmetric sampling as a function of r.
  // Arguments:
  //  etaj : Vector normalized coordinate: etaj = H*posj
  //  etai : Vector normalized coordinate: etai = H*posi
  //  Hdet  : Determinant of the H tensor used to compute eta
  double operator()(const Vector& etaj, const Vector& etai, const Scalar Hdet) const;
  Vector grad(const Vector& etaj, const Vector& etai, const SymTensor& H) const;
  void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;

  // Access our internal data.
  const InterpolatorType& Winterpolator() const;
  const TableKernel<Dim<3>>& baseKernel3d() const;
  const TableKernel<Dim<1>>& baseKernel1d() const;
  Scalar etamax() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  InterpolatorType mInterp;
  TableKernel<Dim<3>> mBaseKernel3d;
  TableKernel<Dim<1>> mBaseKernel1d;  // Only for use with the IdealH algorithm
  Scalar metamax;
};

}

#include "SphericalKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class SphericalKernel;
}

#endif
