//---------------------------------Spheral++----------------------------------//
// DamagedSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of damage.
//
// Created by JMO, Wed Sep  5 14:48:03 PDT 2007
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamagedSoundSpeedPolicy_hh__
#define __Spheral_DamagedSoundSpeedPolicy_hh__

#include "Hydro/SoundSpeedPolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class DamagedSoundSpeedPolicy: public SoundSpeedPolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  DamagedSoundSpeedPolicy();
  virtual ~DamagedSoundSpeedPolicy();
  
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
  DamagedSoundSpeedPolicy(const DamagedSoundSpeedPolicy& rhs);
  DamagedSoundSpeedPolicy& operator=(const DamagedSoundSpeedPolicy& rhs);
};

}

#endif
