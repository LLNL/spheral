//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_ASPHSmooothingScale__
#define __Spheral_NodeSpace_ASPHSmooothingScale__

#include "SmoothingScaleBase.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/CubicHermiteInterpolator.hh"

namespace Spheral {

template<typename Dimension>
class ASPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using InterpolatorType = CubicHermiteInterpolator;

  // Constructors, destructor.
  ASPHSmoothingScale(const TableKernel<Dimension>& W,
                     const Scalar targetNperh,
                     const size_t numPoints = 0u);     // numPoints == 0 ==> use same number of points as TableKernel
  explicit ASPHSmoothingScale();
  ASPHSmoothingScale(const ASPHSmoothingScale& rhs);
  ASPHSmoothingScale& operator=(const ASPHSmoothingScale& rhs);
  virtual ~ASPHSmoothingScale();

  // Time derivative of the smoothing scale.
  virtual
  SymTensor
  smoothingScaleDerivative(const SymTensor& H,
                           const Vector& pos,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh) const override;
  
  // Return a new H, with limiting based on the old value.
  virtual
  SymTensor
  newSmoothingScale(const SymTensor& H,
                    const Vector& pos,
                    const Scalar zerothMoment,
                    const Vector& firstMoment,
                    const SymTensor& secondMomentEta,
                    const SymTensor& secondMomentLab,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    const unsigned nodeListi,
                    const unsigned i) const override;

  // Determine an "ideal" H for the given moments.
  virtual
  SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Vector& pos,
                      const Scalar zerothMoment,
                      const Vector& firstMoment,
                      const SymTensor& secondMomentEta,
                      const SymTensor& secondMomentLab,
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
                      const Scalar nPerh) const override;

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

// We explicitly specialize the time derivatives.
template<> 
Dim<1>::SymTensor
ASPHSmoothingScale<Dim<1> >::smoothingScaleDerivative(const Dim<1>::SymTensor&, 
                                                      const Dim<1>::Vector& pos,
                                                      const Dim<1>::Tensor&,
                                                      const Dim<1>::Scalar hmin,
                                                      const Dim<1>::Scalar hmax,
                                                      const Dim<1>::Scalar hminratio,
                                                      const Dim<1>::Scalar nPerh) const;
template<>
Dim<2>::SymTensor
ASPHSmoothingScale<Dim<2> >::smoothingScaleDerivative(const Dim<2>::SymTensor&, 
                                                      const Dim<2>::Vector& pos,
                                                      const Dim<2>::Tensor&,
                                                      const Dim<2>::Scalar hmin,
                                                      const Dim<2>::Scalar hmax,
                                                      const Dim<2>::Scalar hminratio,
                                                      const Dim<2>::Scalar nPerh) const;
template<>
Dim<3>::SymTensor
ASPHSmoothingScale<Dim<3> >::smoothingScaleDerivative(const Dim<3>::SymTensor&, 
                                                      const Dim<3>::Vector& pos,
                                                      const Dim<3>::Tensor&,
                                                      const Dim<3>::Scalar hmin,
                                                      const Dim<3>::Scalar hmax,
                                                      const Dim<3>::Scalar hminratio,
                                                      const Dim<3>::Scalar nPerh) const;

}

#endif
