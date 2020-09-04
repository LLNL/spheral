//---------------------------------Spheral++----------------------------------//
// EffectiveTensorDamagePolicy -- A ReplaceState policy that does nothing
// when updateAsIncrement is called.
//
// Created by JMO, Thu Dec  4 10:16:18 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_EffectiveTensorDamagePolicy_hh__
#define __Spheral_EffectiveTensorDamagePolicy_hh__

#include "DataBase/ReplaceState.hh"

namespace Spheral {

template<typename Dimension>
class EffectiveTensorDamagePolicy:
    public ReplaceState<Dimension, typename Dimension::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ReplaceState<Dimension, SymTensor>::KeyType KeyType;

  // Constructors, destructor.
  EffectiveTensorDamagePolicy(): ReplaceState<Dimension, SymTensor>() {};
  virtual ~EffectiveTensorDamagePolicy() {};
  
  // Override the updateAsIncrement method to do nothing.
  virtual void updateAsIncrement(const KeyType& /*key*/,
                                 State<Dimension>& /*state*/,
                                 StateDerivatives<Dimension>& /*derivs*/,
                                 const double /*multiplier*/,
                                 const double /*t*/,
                                 const double /*dt*/) {};

private:
  //--------------------------- Private Interface ---------------------------//
  EffectiveTensorDamagePolicy(const EffectiveTensorDamagePolicy& rhs);
  EffectiveTensorDamagePolicy& operator=(const EffectiveTensorDamagePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class EffectiveTensorDamagePolicy;
}

#endif
