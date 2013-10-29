//---------------------------------Spheral++----------------------------------//
// State -- Accumulate and cart the state for a set of physics packages around.
//
// Created by JMO, Thu Aug 26 10:39:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_State_hh__
#define __Spheral_State_hh__

#include <vector>
#include <map>

#include "boost/shared_ptr.hpp"

#include "StateBase.hh"
#include "UpdatePolicyBase.hh"
#include "FieldUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
// template<typename Dimension> class UpdatePolicyBase;
template<typename Dimension> class StateDerivatives;
namespace FieldSpace {
  template<typename Dimension, typename DataType> class Field;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}
namespace PhysicsSpace {
  template<typename Dimension> class Physics;
}

template<typename Dimension>
class State: public StateBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Vector3d Vector3d;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename StateBase<Dimension>::KeyType KeyType;

  typedef std::vector<PhysicsSpace::Physics<Dimension>*> PackageList;
  typedef typename PackageList::iterator PackageIterator;
  typedef typename boost::shared_ptr<UpdatePolicyBase<Dimension> > PolicyPointer;

  // Constructors, destructor.
  State();
  State(DataBaseSpace::DataBase<Dimension>& dataBase,  PackageList& physicsPackages);
  State(DataBaseSpace::DataBase<Dimension>& dataBase,
        PackageIterator physicsPackageBegin,
        PackageIterator physicsPackageEnd);
  State(const State& rhs);
  virtual ~State();

  // Assignment.
  State& operator=(const State& rhs);

  // Override the base method.
  virtual bool operator==(const StateBase<Dimension>& rhs) const;

  // Update the registered state according to the policies.
  void update(StateDerivatives<Dimension>& derivs,
              const double multiplier,
              const double t,
              const double dt);

  // Enroll a policy by itself.
  void enroll(const KeyType& key, PolicyPointer policy);

  // Enroll the given Field and associated update policy
  template<typename Value>
  void enroll(FieldSpace::Field<Dimension, Value>& field, PolicyPointer policy);

  // The base class method for just registering a field.
  virtual void enroll(FieldSpace::FieldBase<Dimension>& field);

  // The base class method for just registering a field list.
  virtual void enroll(FieldSpace::FieldListBase<Dimension>& fieldList);

  // The full set of keys for all policies.
  std::vector<KeyType> policyKeys() const;

  // Return the policy for the specified key.
  PolicyPointer policy(const KeyType& key) const;

  // Return the policy for the specified field.
  template<typename Value>
  PolicyPointer policy(const FieldSpace::Field<Dimension, Value>& field) const;

  // Optionally trip a flag indicating policies should time advance only -- no replacing state!
  // This is useful when you're trying to cheat and reuse derivatives from a prior advance.
  bool timeAdvanceOnly() const;
  void timeAdvanceOnly(const bool x);

private:
  //--------------------------- Private Interface ---------------------------//
  typedef std::map<KeyType, PolicyPointer> PolicyMapType;
  PolicyMapType mPolicyMap;
  bool mTimeAdvanceOnly;
};

}

#include "StateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class State;
}

#endif

