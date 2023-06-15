//---------------------------------Spheral++----------------------------------//
// IncrementSpecificFromTotalPolicy -- policy to update the velocity from the 
//                                  momentum time derivative
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_IncrementSpecificFromTotalPolicy_hh__
#define __Spheral_IncrementSpecificFromTotalPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementSpecificFromTotalPolicy: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:

  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  // Constructors, destructor.

  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0, const std::string& depend1);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~IncrementSpecificFromTotalPolicy();
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:

  const std::string mStateKey;
  const std::string mDerivativeKey;

  //--------------------------- Private Interface ---------------------------//
  IncrementSpecificFromTotalPolicy();
  IncrementSpecificFromTotalPolicy(const IncrementSpecificFromTotalPolicy& rhs);
  IncrementSpecificFromTotalPolicy& operator=(const IncrementSpecificFromTotalPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class IncrementSpecificFromTotalPolicy;
}

#endif
