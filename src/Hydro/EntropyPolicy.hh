//---------------------------------Spheral++----------------------------------//
// EntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the entropy.
//
// Created by JMO, Sun Mar  6 21:57:47 PST 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_EntropyPolicy_hh__
#define __Spheral_EntropyPolicy_hh__

#include <string>

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class Field;
}

template<typename Dimension>
class EntropyPolicy: public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  EntropyPolicy();
  virtual ~EntropyPolicy();
  
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
  EntropyPolicy(const EntropyPolicy& rhs);
  EntropyPolicy& operator=(const EntropyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class EntropyPolicy;
}

#endif
