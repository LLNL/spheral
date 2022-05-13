//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementFieldList -- Update policy for the shearDisplacement pairwise
//                            pairwise field. This one increments a replaced 
//                            field. The replacement reorients the displacment
//                            into the new tangential plane of the particle-
//                            particle contact.
//
// Created by JMP, Sun May 8 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceAndIncrementFieldList_hh__
#define __Spheral_ReplaceAndIncrementFieldList_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename Value>
class ReplaceAndIncrementFieldList: public FieldListUpdatePolicyBase<Dimension, Value> {
public:
  //--------------------------- Public Interface ---------------------------//
  // pull up from the parent
  typedef typename  FieldListUpdatePolicyBase<Dimension, Value>::KeyType KeyType;

  // Constructors, destructor.
  ReplaceAndIncrementFieldList();
  explicit ReplaceAndIncrementFieldList(const std::string& depend0);
  ReplaceAndIncrementFieldList(const std::string& depend0, const std::string& depend1);
  ReplaceAndIncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceAndIncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceAndIncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceAndIncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceAndIncrementFieldList();
  
  static const std::string incrementPrefix() { return "delta "; }
  static const std::string replacePrefix() { return "new "; }

  bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceAndIncrementFieldList(const ReplaceAndIncrementFieldList& rhs);
  ReplaceAndIncrementFieldList& operator=(const ReplaceAndIncrementFieldList& rhs);

};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename Value> class ReplaceAndIncrementFieldList;
}

#endif
