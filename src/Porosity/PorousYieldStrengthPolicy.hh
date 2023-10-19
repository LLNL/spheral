//---------------------------------Spheral++----------------------------------//
// PorousYieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state in the presence of porosity.
//
// Created by JMO, Thu Oct  5 11:04:29 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousYieldStrengthPolicy_hh__
#define __Spheral_PorousYieldStrengthPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class PorousYieldStrengthPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorousYieldStrengthPolicy();
  virtual ~PorousYieldStrengthPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  PorousYieldStrengthPolicy(const PorousYieldStrengthPolicy& rhs);
  PorousYieldStrengthPolicy& operator=(const PorousYieldStrengthPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorousYieldStrengthPolicy;
}

#endif
