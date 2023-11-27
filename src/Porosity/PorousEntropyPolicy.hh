//---------------------------------Spheral++----------------------------------//
// PorousEntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent entropy state accounting for the porosity
// using a P-alpha model
//
// Created by JMO, Thu Oct  5 11:29:39 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousEntropyPolicy_hh__
#define __Spheral_PorousEntropyPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;

template<typename Dimension>
class PorousEntropyPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorousEntropyPolicy();
  virtual ~PorousEntropyPolicy() {}
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  PorousEntropyPolicy(const PorousEntropyPolicy& rhs);
  PorousEntropyPolicy& operator=(const PorousEntropyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorousEntropyPolicy;
}

#endif
