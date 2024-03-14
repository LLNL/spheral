//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScalev2
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Mon Mar 11 10:36:21 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_ASPHSmooothingScalev2__
#define __Spheral_NodeSpace_ASPHSmooothingScalev2__

#include "ASPHSmoothingScale.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/CubicHermiteInterpolator.hh"

namespace Spheral {

template<typename Dimension>
class ASPHSmoothingScalev2: public ASPHSmoothingScale<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using InterpolatorType = CubicHermiteInterpolator;

  // Constructors, destructor.
  ASPHSmoothingScalev2(const TableKernel<Dimension>& W,
                       const Scalar targetNperh,
                       const size_t numPoints = 0u);     // numPoints == 0 ==> use same number of points as TableKernel
  ASPHSmoothingScalev2(const ASPHSmoothingScalev2& rhs);
  ASPHSmoothingScalev2& operator=(const ASPHSmoothingScalev2& rhs);
  virtual ~ASPHSmoothingScalev2();

  // Determine an "ideal" H for the given moments.
  virtual
  SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Vector& pos,
                      const Scalar zerothMoment,
                      const SymTensor& secondMoment,
                      const TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const ConnectivityMap<Dimension>& connectivityMap,
                      const unsigned nodeListi,
                      const unsigned i) const override;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Mesh<Dimension>& mesh,
                      const typename Mesh<Dimension>::Zone& zone,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh) const override { return ASPHSmoothingScale<Dimension>::idealSmoothingScale(H, mesh, zone, hmin, hmax, hminratio, nPerh); }

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values, using the second moment ASPH algorithm
  Scalar equivalentNodesPerSmoothingScale(const Scalar lambdaPsi) const;
  Scalar equivalentLambdaPsi(const Scalar nPerh) const;

  // Access the internal data
  Scalar targetNperh() const                         { return mTargetNperh; }
  Scalar minNperh() const                            { return mMinNperh; }
  Scalar maxNperh() const                            { return mMaxNperh; }
  const InterpolatorType& nPerhInterpolator() const  { return mNperhLookup; }
  const InterpolatorType& WsumInterpolator() const   { return mWsumLookup; }

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mTargetNperh, mMinNperh, mMaxNperh;
  InterpolatorType mNperhLookup, mWsumLookup;
};

}

#endif
