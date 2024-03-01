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
  Scalar kernelValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  Scalar gradValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  Scalar grad2Value(const Scalar etaij, const Scalar Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;
  void kernelAndGradValue(const Scalar etaij, const Scalar Hdet,
                          Scalar& W,
                          Scalar& gW) const;

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const std::vector<Scalar>& etaijs,
                           const std::vector<Scalar>& Hdets,
                           std::vector<Scalar>& kernelValues,
                           std::vector<Scalar>& gradValues) const;

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values.
  Scalar equivalentNodesPerSmoothingScale(const Scalar Wsum) const;

  // Return the equivalent W sum implied by the given number of nodes per smoothing scale.
  Scalar equivalentWsum(const Scalar nPerh) const;

  // Access the internal data
  size_t numPoints() const                               { return mNumPoints; }
  Scalar minNperhLookup() const                          { return mMinNperh; }
  Scalar maxNperhLookup() const                          { return mMaxNperh; }

  // Direct access to our interpolators
  const InterpolatorType& Winterpolator() const          { return mInterp; }
  const InterpolatorType& gradWinterpolator() const      { return mGradInterp; }
  const InterpolatorType& grad2Winterpolator() const     { return mGrad2Interp; }
  const InterpolatorType& nPerhInterpolator() const      { return mNperhLookup; }
  const InterpolatorType& WsumInterpolator() const       { return mWsumLookup; }
  const InterpolatorType& nPerhInterpolatorASPH() const  { return mNperhLookupASPH; }
  const InterpolatorType& WsumInterpolatorASPH() const   { return mWsumLookupASPH; }

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  size_t mNumPoints;
  Scalar mMinNperh, mMaxNperh;
  InterpolatorType mInterp, mGradInterp, mGrad2Interp;  // W, grad W, grad^2 W
  InterpolatorType mNperhLookup, mWsumLookup;           // SPH nperh lookups
  InterpolatorType mNperhLookupASPH, mWsumLookupASPH;   // ASPH nperh lookups
};

}

#include "TableKernelInline.hh"

#endif

