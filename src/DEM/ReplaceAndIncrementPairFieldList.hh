//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementPairFieldList -- Update policy which first replaces the
//                                     pairFieldList in question then increments
//                                     it. Naturally two derivatives fields
//                                     are required one for each of the two
//                                     steps.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceAndIncrementPairFieldList_hh__
#define __Spheral_ReplaceAndIncrementPairFieldList_hh__

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename Value>
class ReplaceAndIncrementPairFieldList: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // pull up from the parent
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  ReplaceAndIncrementPairFieldList(std::initializer_list<std::string> depends = {});
  virtual ~ReplaceAndIncrementPairFieldList() {}
  
  static const std::string incrementPrefix() { return "delta "; }
  static const std::string replacePrefix() { return "new "; }

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
  ReplaceAndIncrementPairFieldList(const ReplaceAndIncrementPairFieldList& rhs);
  ReplaceAndIncrementPairFieldList& operator=(const ReplaceAndIncrementPairFieldList& rhs);

};

}

#endif
