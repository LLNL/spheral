//---------------------------------Spheral++----------------------------------//
// CompositeFieldListPolicy -- An implementation of UpdatePolicyBase which
// consists of a collection of individual Field policies that should match
// the Fields in a FieldList.
//
// Created by JMO, Sun Nov  3 14:11:32 PST 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CompositeFieldListPolicy_hh__
#define __Spheral_CompositeFieldListPolicy_hh__

#include "FieldListUpdatePolicyBase.hh"

#include <vector>
#include <memory> // unique_ptr/shared_ptr

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class CompositeFieldListPolicy: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename std::shared_ptr<UpdatePolicyBase<Dimension> > PolicyPointer;
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  CompositeFieldListPolicy();
  virtual ~CompositeFieldListPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // An alternate method to be called when you want to specify that the derivative information
  // should be assumed to not necessarily be properly time-centered, and therefore you should 
  // only use time advancement ideas, no "replace" or more sophisticated approaches.
  // Default to just calling the generic method.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Add new UpdatePolicies to this thing.
  void push_back(UpdatePolicyBase<Dimension>* policyPtr);

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<std::unique_ptr<UpdatePolicyBase<Dimension>>> mPolicyPtrs;

  CompositeFieldListPolicy(const CompositeFieldListPolicy& rhs);
  CompositeFieldListPolicy& operator=(const CompositeFieldListPolicy& rhs);
};

}

#endif
