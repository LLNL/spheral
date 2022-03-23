//---------------------------------Spheral++----------------------------------//
// ReplaceWithRatioPolicy -- replaces with ratio of two other state fields
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_ReplaceWithRatioPolicy_hh__
#define __Spheral_ReplaceWithRatioPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class ReplaceWithRatioPolicy: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.

  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0, const std::string& depend1);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceWithRatioPolicy(const KeyType& numerator, const KeyType& denomator, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceWithRatioPolicy();
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const KeyType mNumerator;
  const KeyType mDenomenator;

  ReplaceWithRatioPolicy();
  ReplaceWithRatioPolicy(const ReplaceWithRatioPolicy& rhs);
  ReplaceWithRatioPolicy& operator=(const ReplaceWithRatioPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class ReplaceWithRatioPolicy;
}

#endif
