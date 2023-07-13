//---------------------------------Spheral++----------------------------------//
// SphericalRadialKernel
//
// Take a 1D Kernel and build a specialized 1D tabulated version appropriate
// for use with the radial spherical SPH algorithm
//
// Created by JMO, Wed Jul  5 10:46:14 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalRadialKernel_hh__
#define __Spheral_SphericalRadialKernel_hh__

#include "TableKernel.hh"
#include "Utilities/CubicHermiteInterpolator.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalRadialKernel {

public:
  //--------------------------- Public Interface ---------------------------//
  using InterpolatorType = CubicHermiteInterpolator;
  using Scalar = Dim<1>::Scalar;
  using Vector = Dim<1>::Vector;
  using Tensor = Dim<1>::Tensor;
  using SymTensor = Dim<1>::SymTensor;

  // Constructor.
  // Takes a normal 1D TableKernel and constructs the integral form appropriate
  // for 1D spherical coordinates.
  template<typename KernelType>
  explicit
  SphericalRadialKernel(const KernelType& kernel,
                        const unsigned numIntegral = 5000u,
                        const unsigned numKernel = 200u,
                        const bool useInterpolation = true);
  SphericalRadialKernel(const SphericalRadialKernel& rhs);

  // Destructor.
  virtual ~SphericalRadialKernel();

  // Assignment.
  SphericalRadialKernel& operator=(const SphericalRadialKernel& rhs);

  // Comparisons
  bool operator==(const SphericalRadialKernel& rhs) const;

  // These methods taking a Vector eta and Vector position are the special methods
  // allowing this kernel to implement the asymmetric sampling as a function of r.
  // Arguments:
  //  etaj : Vector normalized coordinate: etaj = H*posj
  //  etai : Vector normalized coordinate: etai = H*posi
  //  Hdet  : Determinant of the H tensor used to compute eta
  Scalar operator()(const Vector& etaj, const Vector& etai, const Scalar Hdet) const;
  Vector grad(const Vector& etaj, const Vector& etai, const SymTensor& H) const;
  void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;

  // Look up/compute the volume nornalization (using interpolation or not based on switch)
  Scalar volumeNormalization(const Scalar eta) const;
  Scalar gradAInv(const Scalar eta) const;

  // Access our internal data.
  const InterpolatorType& Ainterpolator() const;
  const InterpolatorType& gradAInvInterpolator() const;
  const TableKernel<Dim<1>>& baseKernel1d() const;
  Scalar etamax() const;
  Scalar etacutoff() const;
  bool useInterpolation() const;
  void useInterpolation(const bool x);

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  InterpolatorType mAInterp, mGradAInvInterp;
  TableKernel<Dim<1>> mBaseKernel;
  Scalar metamax, metacutoff, mAetacutoff, mGradAInvEtacutoff;
  unsigned mNumIntegral;
  bool mUseInterpolation;
};

}

#include "SphericalRadialKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class SphericalRadialKernel;
}

#endif

