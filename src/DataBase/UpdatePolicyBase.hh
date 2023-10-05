//---------------------------------Spheral++----------------------------------//
// UpdatePolicyBase -- Base/interface class for the policies on how state 
// variables are to be updated.
//
// Note that we currenly only support specifying [0,6] dependencies.  There's 
// no problem adding more constructors to expand the allowed number of 
// dependencies, this is just where I've arbitrarily drawn the line for now.
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_UpdatePolicyBase_hh__
#define __Spheral_UpdatePolicyBase_hh__

#include <string>
#include <vector>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class UpdatePolicyBase {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef std::string KeyType;

  // Constructors, destructor.
  UpdatePolicyBase();
  explicit UpdatePolicyBase(const std::string& depend0);
  UpdatePolicyBase(const std::string& depend0, 
                   const std::string& depend1);
  UpdatePolicyBase(const std::string& depend0, 
                   const std::string& depend1,
                   const std::string& depend2);
  UpdatePolicyBase(const std::string& depend0, 
                   const std::string& depend1, 
                   const std::string& depend2, 
                   const std::string& depend3);
  UpdatePolicyBase(const std::string& depend0, 
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3, 
                   const std::string& depend4);
  UpdatePolicyBase(const std::string& depend0, 
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4, 
                   const std::string& depend5);
  UpdatePolicyBase(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2, 
                   const std::string& depend3, 
                   const std::string& depend4,
                   const std::string& depend5,
                   const std::string& depend6);
  virtual ~UpdatePolicyBase() {};
  
  // The methods you overload to define behavior mapping new state field 
  // from derivative methods.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) = 0;

  // An alternate method to be called when you want to specify that the derivative information
  // should be assumed to not necessarily be properly time-centered, and therefore you should 
  // only use time advancement ideas, no "replace" or more sophisticated approaches.
  // Default to just calling the generic method.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    this->update(key, state, derivs, multiplier, t, dt);
  }

  // Require descendents to define an equivalence operator.
  virtual bool operator==(const UpdatePolicyBase& rhs) const = 0;
  bool operator!=(const UpdatePolicyBase& rhs) const;

  // Should this policy be cloned per Field when registering for a FieldList?
  virtual bool clonePerField() const { return false; }

  // Test is this policy is for independent or dependent state.
  bool independent() const;
  bool dependent() const;

  // Return the set of field names that this state depends upon (if any).
  const std::vector<std::string>& dependencies() const;

  // Allow the addition of new dependencies.
  void addDependency(const std::string& depend);

  // The wildcard string for comparing dependency keys.
  static const std::string wildcard() { return "*"; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<std::string> mDependencies;

  UpdatePolicyBase(const UpdatePolicyBase& rhs);
  UpdatePolicyBase& operator=(const UpdatePolicyBase& rhs);
};

}

#ifndef __GCCXML__
#include "UpdatePolicyBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class UpdatePolicyBase;
}

#endif
