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

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename Value>
class ReplaceAndIncrementPairFieldList: public FieldListUpdatePolicyBase<Dimension, Value> {
public:
  //--------------------------- Public Interface ---------------------------//
  // pull up from the parent
  typedef typename  FieldListUpdatePolicyBase<Dimension, Value>::KeyType KeyType;

  // Constructors, destructor.
  ReplaceAndIncrementPairFieldList();
  explicit ReplaceAndIncrementPairFieldList(const std::string& depend0);
  ReplaceAndIncrementPairFieldList(const std::string& depend0, const std::string& depend1);
  ReplaceAndIncrementPairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceAndIncrementPairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceAndIncrementPairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceAndIncrementPairFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceAndIncrementPairFieldList();
  
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
