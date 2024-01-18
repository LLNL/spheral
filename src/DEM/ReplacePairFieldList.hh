//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementPairFieldList -- Update policy which replaces the values
//                                     of a pairFieldList.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplacePairFieldList_hh__
#define __Spheral_ReplacePairFieldList_hh__

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class ReplacePairFieldList: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // pull up from the parent
  using KeyType = typename  UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  ReplacePairFieldList(std::initializer_list<std::string> depends = {});
  virtual ~ReplacePairFieldList() {}
  
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

#endif
