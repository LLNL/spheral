//---------------------------------Spheral++----------------------------------//
// GammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_GammaPolicy_hh__
#define __Spheral_GammaPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class GammaPolicy: public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  GammaPolicy();
  virtual ~GammaPolicy();
  
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
  GammaPolicy(const GammaPolicy& rhs);
  GammaPolicy& operator=(const GammaPolicy& rhs);
};

}

#endif
