//---------------------------------Spheral++----------------------------------//
// ReplacePairFieldList -- An implementation of FieldListUpdatePolicyBase appropriate for
// when 'ya just want to Replace by derivatives:  x1 = x0 + A*dx/dt
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplacePairFieldList_hh__
#define __Spheral_ReplacePairFieldList_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class ReplacePairFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // pull up from the parent
  typedef typename  FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  ReplacePairFieldList();
  explicit ReplacePairFieldList(const std::string& depend0);
  ReplacePairFieldList(const std::string& depend0, const std::string& depend1);
  ReplacePairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplacePairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplacePairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplacePairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplacePairFieldList();
  
  static const std::string prefix() { return "new "; }
  
  bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

private:
  //--------------------------- Private Interface ---------------------------//
  ReplacePairFieldList(const ReplacePairFieldList& rhs);
  ReplacePairFieldList& operator=(const ReplacePairFieldList& rhs);

};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class ReplacePairFieldList;
}

#endif
