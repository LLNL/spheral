//---------------------------------Spheral++----------------------------------//
// FieldUpdatePolicy -- Base/interface class for the policies on how 
// Field state variables are to be updated.
//
// Created by JMO, Sun Feb 13 20:50:53 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldUpdatePolicy_hh__
#define __Spheral_FieldUpdatePolicy_hh__

#include "UpdatePolicyBase.hh"

#include <memory> // unique_ptr/shared_ptr

namespace Spheral {

template<typename Dimension>
class FieldUpdatePolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  FieldUpdatePolicy(std::initializer_list<std::string> depends = {});
  virtual ~FieldUpdatePolicy() {};
  
  // Should this policy be cloned per Field when registering for a FieldList?
  // Setting this to true is the only purpose of this class, so that policies
  // intended for Fields can be easily registered on FieldLists.  Thie results
  // in each Field in the FieldList being registered separately with copies of
  // the shared_ptr to the policy.
  virtual bool clonePerField() const override { return true; }

private:
  //--------------------------- Private Interface ---------------------------//
  FieldUpdatePolicy(const FieldUpdatePolicy& rhs);
  FieldUpdatePolicy& operator=(const FieldUpdatePolicy& rhs);
};

}

#include "FieldUpdatePolicyInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class FieldUpdatePolicy;
}

#endif
