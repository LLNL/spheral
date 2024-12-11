//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_FiniteVolumeViscosity__
#define __Spheral_FiniteVolumeViscosity__

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension, typename Value> class PairwiseField;

template<typename Dimension>
class FiniteVolumeViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ConstBoundaryIterator = typename ArtificialViscosity<Dimension>::ConstBoundaryIterator;
  using PairQPiType = PairwiseField<Dimension, std::pair<Scalar, Scalar>>;

  // Constructor, destructor
  FiniteVolumeViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const TableKernel<Dimension>& WT);
  virtual ~FiniteVolumeViscosity() = default;

  // No default construction, copying, or assignment
  FiniteVolumeViscosity() = delete;
  FiniteVolumeViscosity(const FiniteVolumeViscosity&) = delete;
  FiniteVolumeViscosity& operator=(const FiniteVolumeViscosity&) const = delete;

  // We compute a scalar QPi
  virtual std::type_index QPiType()                  const override { return std::type_index(typeid(Scalar)); }

  // We are going to use a velocity gradient
  virtual bool requireVelocityGradient()             const override { return true; }

  // Add our contribution to the derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  // Override the method of computing the velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs) override;

  // Restart methods.
  virtual std::string label()                        const override { return "FiniteVolumeViscosity"; }
};

}

#endif
