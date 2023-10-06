//---------------------------------Spheral++----------------------------------//
// PureReplaceState -- replaces one field values with those of 
//                     another specified by its key
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_PureReplaceState_hh__
#define __Spheral_PureReplaceState_hh__

#include "DataBase/ReplaceState.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class PureReplaceState: public ReplaceState<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename ReplaceState<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  explicit PureReplaceState(const KeyType& derivFieldListKey);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~PureReplaceState();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Override the base method for computing the new field key.
  // This is the only tweak this policy imposes over the standard ReplaceState.
  virtual KeyType replaceStateKey(const KeyType& fkey) const override { return mReplaceKey; }

private:
  //--------------------------- Private Interface ---------------------------//
  const KeyType mReplaceKey;

  PureReplaceState();
  PureReplaceState(const PureReplaceState& rhs);
  PureReplaceState& operator=(const PureReplaceState& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class PureReplaceState;
}

#endif
