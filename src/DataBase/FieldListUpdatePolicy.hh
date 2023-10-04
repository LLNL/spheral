//---------------------------------Spheral++----------------------------------//
// FieldListUpdatePolicyBase -- Base/interface class for the policies on how 
// FieldList state variables are to be updated.
//
// Created by JMO, Sun Feb 13 20:50:53 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListUpdatePolicyBase_hh__
#define __Spheral_FieldListUpdatePolicyBase_hh__

#include "UpdatePolicyBase.hh"

#include <string>
#include <map>
#include <memory> // unique_ptr/shared_ptr

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldListUpdatePolicyBase: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  FieldListUpdatePolicyBase();
  explicit FieldListUpdatePolicyBase(const std::string& depend0);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1, 
                            const std::string& depend2, 
                            const std::string& depend3);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2,
                            const std::string& depend3, 
                            const std::string& depend4);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2,
                            const std::string& depend3,
                            const std::string& depend4, 
                            const std::string& depend5);
  FieldListUpdatePolicyBase(const std::string& depend0,
                            const std::string& depend1,
                            const std::string& depend2, 
                            const std::string& depend3, 
                            const std::string& depend4,
                            const std::string& depend5,
                            const std::string& depend6);
  virtual ~FieldListUpdatePolicyBase() {};
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // An alternate method to be called when you want to specify that the derivative information
  // should be assumed to not necessarily be properly time-centered, and therefore you should 
  // only use time advancement ideas, no "replace" or more sophisticated approaches.
  // Default to just calling the generic method.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Set the policy for a given NodeList
  void enroll(const std::string& nodeListName, std::shared_ptr<UpdatePolicyBase<Dimension>> policy);
  bool havePolicyForNodeList(const std::string& nodeListName) const;
  const std::shared_ptr<UpdatePolicyBase<Dimension>>& policyForNodeList(const std::string& nodeListName) const;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldListUpdatePolicyBase(const FieldListUpdatePolicyBase& rhs);
  FieldListUpdatePolicyBase& operator=(const FieldListUpdatePolicyBase& rhs);

  // The set of UpdatePolicies by NodeList.
  std::map<std::string, std::shared_ptr<UpdatePolicyBase<Dimension>>> mNodeListPolicies;
};

}

#include "FieldListUpdatePolicyBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename DataType> class FieldListUpdatePolicyBase;
}

#endif
