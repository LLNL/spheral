//---------------------------------Spheral++----------------------------------//
// SoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent sound speed state.
//
// Created by JMO, Wed Sep 15 10:47:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_SoundSpeedPolicy_hh__
#define __Spheral_SoundSpeedPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class SoundSpeedPolicy: public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  SoundSpeedPolicy();
  virtual ~SoundSpeedPolicy();
  
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
  SoundSpeedPolicy(const SoundSpeedPolicy& rhs);
  SoundSpeedPolicy& operator=(const SoundSpeedPolicy& rhs);
};

}

#endif
