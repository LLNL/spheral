//---------------------------------Spheral++----------------------------------//
// PressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent pressure state.
//
// Created by JMO, Wed Sep 15 10:47:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PressurePolicy_hh__
#define __Spheral_PressurePolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class PressurePolicy: public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  PressurePolicy();
  virtual ~PressurePolicy();
  
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
  PressurePolicy(const PressurePolicy& rhs);
  PressurePolicy& operator=(const PressurePolicy& rhs);
};

}

#endif
