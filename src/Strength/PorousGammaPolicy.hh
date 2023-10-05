//---------------------------------Spheral++----------------------------------//
// PorousGammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state in the presence of porosity.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousGammaPolicy_hh__
#define __Spheral_PorousGammaPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class PorousGammaPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorousGammaPolicy();
  virtual ~PorousGammaPolicy();
  
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
  PorousGammaPolicy(const PorousGammaPolicy& rhs);
  PorousGammaPolicy& operator=(const PorousGammaPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorousGammaPolicy;
}

#endif
