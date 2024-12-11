//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosity__
#define __Spheral_MonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension, typename Value> class PairwiseField;

template<typename Dimension>
class MonaghanGingoldViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ConstBoundaryIterator = typename ArtificialViscosity<Dimension>::ConstBoundaryIterator;
  using PairQPiType = PairwiseField<Dimension, std::pair<Scalar, Scalar>>;

  // Constructors.
  MonaghanGingoldViscosity(const Scalar Clinear,
                           const Scalar Cquadratic,
                           const bool linearInExpansion,
                           const bool quadraticInExpansion,
                           const TableKernel<Dimension>& kernel);
  virtual ~MonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  MonaghanGingoldViscosity() = delete;
  MonaghanGingoldViscosity(const MonaghanGingoldViscosity&) = delete;
  MonaghanGingoldViscosity& operator=(const MonaghanGingoldViscosity&) = delete;

  // We compute a Scalar for QPi
  virtual std::type_index QPiType()      const override { return std::type_index(typeid(Scalar)); }

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Add our contribution to the derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  // Restart methods.
  virtual std::string label()             const { return "MonaghanGingoldViscosity"; }

  // Access data members
  bool linearInExpansion()                const { return mLinearInExpansion; }
  bool quadraticInExpansion()             const { return mQuadraticInExpansion; }
  void linearInExpansion(const bool x)          { mLinearInExpansion = x; }
  void quadraticInExpansion(const bool x)       { mQuadraticInExpansion = x; }

  const PairQPiType& pairQPi()            const { VERIFY2(mPairQPiPtr, "MonaghanGingoldViscosity pairQPi unintialized");  return *mPairQPiPtr; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion, mQuadraticInExpansion;
  std::unique_ptr<PairQPiType> mPairQPiPtr;
};

}

#endif
