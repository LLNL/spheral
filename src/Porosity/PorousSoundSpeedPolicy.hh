//---------------------------------Spheral++----------------------------------//
// PorousSoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent sound speed state accounting for the porosity
// using a P-alpha model
//
// Created by JMO, Thu Sep 28 16:39:01 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousSoundSpeedPolicy_hh__
#define __Spheral_PorousSoundSpeedPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;

template<typename Dimension>
class PorousSoundSpeedPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorousSoundSpeedPolicy();
  virtual ~PorousSoundSpeedPolicy();
  
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
  PorousSoundSpeedPolicy(const PorousSoundSpeedPolicy& rhs);
  PorousSoundSpeedPolicy& operator=(const PorousSoundSpeedPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorousSoundSpeedPolicy;
}

#endif
