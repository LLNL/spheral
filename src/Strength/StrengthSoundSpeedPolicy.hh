//---------------------------------Spheral++----------------------------------//
// StrengthSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of strength.
//
// Created by JMO, Thu Feb  2 14:08:37 PST 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrengthSoundSpeedPolicy_hh__
#define __Spheral_StrengthSoundSpeedPolicy_hh__

#include "Hydro/SoundSpeedPolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class StrengthSoundSpeedPolicy: public SoundSpeedPolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  StrengthSoundSpeedPolicy();
  virtual ~StrengthSoundSpeedPolicy();
  
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
  StrengthSoundSpeedPolicy(const StrengthSoundSpeedPolicy& rhs);
  StrengthSoundSpeedPolicy& operator=(const StrengthSoundSpeedPolicy& rhs);
};

}

#endif
