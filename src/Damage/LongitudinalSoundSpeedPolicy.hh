//---------------------------------Spheral++----------------------------------//
// LongitudinalSoundSpeedPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the dependent longitudinal sound speed.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_LongitudinalSoundSpeedPolicy_hh__
#define __Spheral_LongitudinalSoundSpeedPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class LongitudinalSoundSpeedPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  LongitudinalSoundSpeedPolicy();
  virtual ~LongitudinalSoundSpeedPolicy();
  
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
  LongitudinalSoundSpeedPolicy(const LongitudinalSoundSpeedPolicy& rhs);
  LongitudinalSoundSpeedPolicy& operator=(const LongitudinalSoundSpeedPolicy& rhs);
};

}

#endif
