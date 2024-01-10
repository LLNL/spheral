//---------------------------------Spheral++----------------------------------//
// IncrementSpecificFromTotalPolicy -- policy to update the velocity from the 
//                                  momentum time derivative
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_IncrementSpecificFromTotalPolicy_hh__
#define __Spheral_IncrementSpecificFromTotalPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementSpecificFromTotalPolicy: public UpdatePolicyBase<Dimension> {
public:

  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  IncrementSpecificFromTotalPolicy(std::initializer_list<std::string> depends, const std::string& stateKey, const std::string& derivKey);
  IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey);
   ~IncrementSpecificFromTotalPolicy();
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }
  
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
