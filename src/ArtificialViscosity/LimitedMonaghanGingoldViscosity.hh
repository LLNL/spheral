//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimitedMonaghanGingoldViscosity__
#define __Spheral_LimitedMonaghanGingoldViscosity__

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class LimitedMonaghanGingoldViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using ConstBoundaryIterator = typename ArtificialViscosity<Dimension>::ConstBoundaryIterator;
  using PairQPiType = PairwiseField<Dimension, std::pair<Scalar, Scalar>>;

  // Constructors.
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac,
                                  const TableKernel<Dimension>& kernel);
  virtual ~LimitedMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  LimitedMonaghanGingoldViscosity() = delete;
  LimitedMonaghanGingoldViscosity(const LimitedMonaghanGingoldViscosity&) = delete;
  LimitedMonaghanGingoldViscosity& operator=(const LimitedMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Add our contribution to the derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;


  // Access our data
  Scalar etaCritFrac()              const { return mEtaCritFrac; }
  Scalar etaFoldFrac()              const { return mEtaFoldFrac; }

  void etaCritFrac(const Scalar x)        { mEtaCritFrac = x; }
  void etaFoldFrac(const Scalar x)        { mEtaFoldFrac = x; }

  // Restart methods.
  virtual std::string label()       const { return "LimitedMonaghanGingoldViscosity"; }

protected:
  //--------------------------- Private Interface ---------------------------//
  double mEtaCritFrac, mEtaFoldFrac;

  using MonaghanGingoldViscosity<Dimension>::mLinearInExpansion;
  using MonaghanGingoldViscosity<Dimension>::mQuadraticInExpansion;
};

}

#endif
