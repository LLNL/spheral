//---------------------------------Spheral++----------------------------------//
// CompositeFieldListPolicy -- An implementation of UpdatePolicyBase which
// consists of a collection of individual Field policies that should match
// the Fields in a FieldList.
//
// Created by JMO, Sun Nov  3 14:11:32 PST 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CompositeFieldListPolicy_hh__
#define __Spheral_CompositeFieldListPolicy_hh__

#include <vector>
#include "boost/shared_ptr.hpp"
#include "boost/ptr_container/ptr_vector.hpp"
#include "FieldListUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class CompositeFieldListPolicy: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename boost::shared_ptr<UpdatePolicyBase<Dimension> > PolicyPointer;
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

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Add new UpdatePolicies to this thing.
  void push_back(UpdatePolicyBase<Dimension>* policyPtr);

private:
  //--------------------------- Private Interface ---------------------------//
  boost::ptr_vector<UpdatePolicyBase<Dimension> > mPolicyPtrs;

  CompositeFieldListPolicy(const CompositeFieldListPolicy& rhs);
  CompositeFieldListPolicy& operator=(const CompositeFieldListPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class CompositeFieldListPolicy;
}

#endif
