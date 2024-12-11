//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//
// Created by J. Michael Owen, Mon Sep  2 14:45:35 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorMonaghanGingoldViscosity__
#define __Spheral_TensorMonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"

#include <memory>

namespace Spheral {

template<typename Dimension, typename Value> class PairwiseField;

template<typename Dimension>
class TensorMonaghanGingoldViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using PairQPiType = PairwiseField<Dimension, std::pair<Tensor, Tensor>>;

  // Constructors and destuctor
  TensorMonaghanGingoldViscosity(const Scalar Clinear,
                                 const Scalar Cquadratic,
                                 const TableKernel<Dimension>& kernel);
  ~TensorMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  TensorMonaghanGingoldViscosity() = delete;
  TensorMonaghanGingoldViscosity(const TensorMonaghanGingoldViscosity&) = delete;
  TensorMonaghanGingoldViscosity& operator=(const TensorMonaghanGingoldViscosity&) = delete;

  // We compute a Tensor for QPi
  virtual std::type_index QPiType()      const override { return std::type_index(typeid(Tensor)); }

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

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
  virtual std::string label() const override { return "TensorMonaghanGingoldViscosity"; }

  // Access data members
  const PairQPiType& pairQPi()         const { VERIFY2(mPairQPiPtr, "TensorMonaghanGingoldViscosity pairQPi unintialized");  return *mPairQPiPtr; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairQPiType> mPairQPiPtr;
};

}

#endif
