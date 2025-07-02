//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#include "Kernel.hh"
#include "Utilities/QuadraticInterpolator.hh"
#include "Utilities/CubicHermiteInterpolator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class TableKernel: public Kernel<Dimension, TableKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using InterpolatorType = QuadraticInterpolator;
  using NperhInterpolatorType = CubicHermiteInterpolator;

  // Constructors.
  template<typename KernelType>
  TableKernel(const KernelType& kernel,
              const unsigned numPoints = 100u,
              const Scalar minNperh = 0.25,
              const Scalar maxNperh = 64.0);
  TableKernel(const TableKernel<Dimension>& rhs);

  // Destructor.
  virtual ~TableKernel();

  // Assignment.
  TableKernel& operator=(const TableKernel& rhs);

  // Equivalence
  bool operator==(const TableKernel& rhs) const;

  // Return the kernel weight for a given normalized distance or position.
  RAJA_HOST_DEVICE Scalar kernelValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  RAJA_HOST_DEVICE Scalar gradValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  RAJA_HOST_DEVICE Scalar grad2Value(const Scalar etaij, const Scalar Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  RAJA_HOST_DEVICE void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;
  RAJA_HOST_DEVICE void kernelAndGradValue(const Scalar etaij, const Scalar Hdet,
                          Scalar& W,
                          Scalar& gW) const;

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const Scalar* etaijs,
                           const Scalar* Hdets,
                           Scalar* kernelValues,
                           Scalar* gradValues,
                           const size_t n) const;

  // Special kernel values for use in finding smoothing scales (SPH and ASPH versions)
  // ***These are only intended for use adapting smoothing scales***, and are used
  // for the succeeding equivalentNodesPerSmoothingScale lookups!
  RAJA_HOST_DEVICE Scalar kernelValueSPH(const Scalar etaij) const;
  RAJA_HOST_DEVICE Scalar kernelValueASPH(const Scalar etaij, const Scalar nPerh) const;

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values, using the zeroth moment SPH algorithm
  RAJA_HOST_DEVICE Scalar equivalentNodesPerSmoothingScale(const Scalar Wsum) const;
  RAJA_HOST_DEVICE Scalar equivalentWsum(const Scalar nPerh) const;

  // Access the internal data
  RAJA_HOST_DEVICE size_t numPoints() const                                    { return mNumPoints; }
  RAJA_HOST_DEVICE Scalar minNperhLookup() const                               { return mMinNperh; }
  RAJA_HOST_DEVICE Scalar maxNperhLookup() const                               { return mMaxNperh; }

  // Direct access to our interpolators
  RAJA_HOST_DEVICE const InterpolatorType& Winterpolator() const               { return mInterp; }
  RAJA_HOST_DEVICE const InterpolatorType& gradWinterpolator() const           { return mGradInterp; }
  RAJA_HOST_DEVICE const InterpolatorType& grad2Winterpolator() const          { return mGrad2Interp; }
  RAJA_HOST_DEVICE const NperhInterpolatorType& nPerhInterpolator() const      { return mNperhLookup; }
  RAJA_HOST_DEVICE const NperhInterpolatorType& WsumInterpolator() const       { return mWsumLookup; }

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  size_t mNumPoints;
  Scalar mTargetNperh, mMinNperh, mMaxNperh;
  InterpolatorType mInterp, mGradInterp, mGrad2Interp;       // W, grad W, grad^2 W
  NperhInterpolatorType mNperhLookup, mWsumLookup;           // SPH nperh lookups

  using Kernel<Dimension, TableKernel<Dimension>>::mVolumeNormalization;
  using Kernel<Dimension, TableKernel<Dimension>>::mKernelExtent;
  using Kernel<Dimension, TableKernel<Dimension>>::mInflectionPoint;
};

}

#include "TableKernelInline.hh"

#endif

