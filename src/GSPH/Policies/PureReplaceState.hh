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
  PureReplaceState(const KeyType& derivFieldListKey,
                   std::initializer_list<std::string> depends = {});
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

  static const std::string prefix() { return "new "; }

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
