//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_YieldStrengthPolicy_hh__
#define __Spheral_YieldStrengthPolicy_hh__

#include <string>

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class YieldStrengthPolicy: 
    public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  YieldStrengthPolicy();
  virtual ~YieldStrengthPolicy();
  
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
  YieldStrengthPolicy(const YieldStrengthPolicy& rhs);
  YieldStrengthPolicy& operator=(const YieldStrengthPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class YieldStrengthPolicy;
}

#endif
