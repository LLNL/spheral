//---------------------------------Spheral++----------------------------------//
// PureReplaceFieldList -- replaces one fieldlists values with those of 
//                         another specified by its key
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_PureReplaceFieldList_hh__
#define __Spheral_PureReplaceFieldList_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class PureReplaceFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.

  explicit PureReplaceFieldList(const KeyType& derivFieldListKey);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  PureReplaceFieldList(const KeyType& derivFieldListKey, const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~PureReplaceFieldList();
  
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
  const KeyType mReplaceKey;

  PureReplaceFieldList();
  PureReplaceFieldList(const PureReplaceFieldList& rhs);
  PureReplaceFieldList& operator=(const PureReplaceFieldList& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class PureReplaceFieldList;
}

#endif
