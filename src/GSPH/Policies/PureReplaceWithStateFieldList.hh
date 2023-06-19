//---------------------------------Spheral++----------------------------------//
// PureReplaceWithStateFieldList -- replaces one fieldlists values with  
//                                  a state field list specified by its key
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_PureReplaceWithStateFieldList_hh__
#define __Spheral_PureReplaceWithStateFieldList_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class PureReplaceWithStateFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.

  explicit PureReplaceWithStateFieldList(const KeyType& derivFieldListKey);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  PureReplaceWithStateFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~PureReplaceWithStateFieldList();
  
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

  PureReplaceWithStateFieldList();
  PureReplaceWithStateFieldList(const PureReplaceWithStateFieldList& rhs);
  PureReplaceWithStateFieldList& operator=(const PureReplaceWithStateFieldList& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class PureReplaceWithStateFieldList;
}

#endif
