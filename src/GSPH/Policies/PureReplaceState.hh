//---------------------------------Spheral++----------------------------------//
// PureReplaceState -- replaces one field values with those of 
//                     another specified by its key
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_PureReplaceState_hh__
#define __Spheral_PureReplaceState_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class PureReplaceState: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  explicit PureReplaceState(const KeyType& derivFieldListKey);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  PureReplaceState(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~PureReplaceState();
  
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
