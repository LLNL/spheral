//---------------------------------Spheral++----------------------------------//
// ReplaceWithRatioPolicy -- replaces with ratio of two other state fields
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_ReplaceWithRatioPolicy_hh__
#define __Spheral_ReplaceWithRatioPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class ReplaceWithRatioPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator);
  ReplaceWithRatioPolicy(std::initializer_list<std::string> depends,
                         const KeyType& numerator, const KeyType& denomator);
  virtual ~ReplaceWithRatioPolicy();
  
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
  //--------------------------- Private Interface ---------------------------//
  const KeyType mNumerator;
  const KeyType mDenomenator;

  ReplaceWithRatioPolicy();
  ReplaceWithRatioPolicy(const ReplaceWithRatioPolicy& rhs);
  ReplaceWithRatioPolicy& operator=(const ReplaceWithRatioPolicy& rhs);
};

}

#endif
