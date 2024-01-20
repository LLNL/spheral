//---------------------------------Spheral++----------------------------------//
// DamageGradientPolicy -- A ReplaceState policy that does nothing
// when updateAsIncrement is called.
//
// Created by JMO, Thu Dec  4 10:16:18 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamageGradientPolicy_hh__
#define __Spheral_DamageGradientPolicy_hh__

#include "DataBase/ReplaceState.hh"

namespace Spheral {

template<typename Dimension>
class DamageGradientPolicy:
    public ReplaceState<Dimension, typename Dimension::Vector> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Vector Vector;
  typedef typename ReplaceState<Dimension, Vector>::KeyType KeyType;

  // Constructors, destructor.
  DamageGradientPolicy(): ReplaceState<Dimension, Vector>() {};
  virtual ~DamageGradientPolicy() {};
  
  // Override the updateAsIncrement method to do nothing.
  virtual void updateAsIncrement(const KeyType& /*key*/,
                                 State<Dimension>& /*state*/,
                                 StateDerivatives<Dimension>& /*derivs*/,
                                 const double /*multiplier*/,
                                 const double /*t*/,
                                 const double /*dt*/) {};

private:
  //--------------------------- Private Interface ---------------------------//
  DamageGradientPolicy(const DamageGradientPolicy& rhs);
  DamageGradientPolicy& operator=(const DamageGradientPolicy& rhs);
};

}

#endif
