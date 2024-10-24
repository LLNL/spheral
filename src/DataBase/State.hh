//---------------------------------Spheral++----------------------------------//
// State -- Accumulate and cart the state for a set of physics packages around.
//
// Created by JMO, Thu Aug 26 10:39:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_State_hh__
#define __Spheral_State_hh__

#include "StateBase.hh"
#include "UpdatePolicyBase.hh"

#include <vector>
#include <map>
#include <memory>

namespace Spheral {

// Forward declarations.
// template<typename Dimension> class UpdatePolicyBase;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class DataBase;
template<typename Dimension> class Physics;

template<typename Dimension>
class State: public StateBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Vector3d = typename Dimension::Vector3d;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using KeyType = typename StateBase<Dimension>::KeyType;

  using PackageList = std::vector<Physics<Dimension>*>;
  using PackageIterator = typename PackageList::iterator;
  using PolicyPointer = typename std::shared_ptr<UpdatePolicyBase<Dimension>>;

  // Constructors, destructor.
  State();
  State(DataBase<Dimension>& dataBase,  PackageList& physicsPackages);
  State(DataBase<Dimension>& dataBase,
        PackageIterator physicsPackageBegin,
        PackageIterator physicsPackageEnd);
  State(const State& rhs);
  virtual ~State();

  // Assignment.
  State& operator=(const State& rhs);

  // Override the base method.
  virtual bool operator==(const StateBase<Dimension>& rhs) const override;

  // Update the registered state according to the policies.
  void update(StateDerivatives<Dimension>& derivs,
              const double multiplier,
              const double t,
              const double dt);

  // Enroll a policy by itself.
  void enroll(const KeyType& key, PolicyPointer policy);

  // Remove a policy.
  void removePolicy(const KeyType& key);
  void removePolicy(FieldBase<Dimension>& field);
  void removePolicy(FieldListBase<Dimension>& field,
                    const bool clonePerField);

  // Enroll the given Field and associated update policy
  void enroll(FieldBase<Dimension>& field, PolicyPointer policy);

  // Enroll the given FieldList and associated update policy
  // This method queries the "clonePerField" method of the policy, and
  // if true enrolls each Field in the FieldList with a copy of the policy.
  // Otherwise the FieldList is enrolled directly as normal, and the policy is
  // assumed to handle a FieldList directly.
  void enroll(FieldListBase<Dimension>& fieldList, PolicyPointer policy);

  // The base class method for just registering a field.
  virtual void enroll(FieldBase<Dimension>& field) override;
  virtual void enroll(std::shared_ptr<FieldBase<Dimension>>& fieldPtr) override;

  // The base class method for just registering a field list.
  virtual void enroll(FieldListBase<Dimension>& fieldList) override;

  // The full set of keys for all policies.
  std::vector<KeyType> policyKeys() const;

  // Return the policy for the specified key.
  PolicyPointer policy(const KeyType& key) const;

  // Return all the policies associated with the given Field Key
  std::map<KeyType, PolicyPointer> policies(const KeyType& fieldKey) const;

  // Return the policy for the specified field.
  template<typename Value>
  PolicyPointer policy(const Field<Dimension, Value>& field) const;

  // Optionally trip a flag indicating policies should time advance only -- no replacing state!
  // This is useful when you're trying to cheat and reuse derivatives from a prior advance.
  bool timeAdvanceOnly() const;
  void timeAdvanceOnly(const bool x);

private:
  //--------------------------- Private Interface ---------------------------//
  using PolicyMapType = std::map<KeyType, std::map<KeyType, PolicyPointer>>;
  PolicyMapType mPolicyMap;
  bool mTimeAdvanceOnly;
};

}

#include "StateInline.hh"

#endif

