//---------------------------------Spheral++----------------------------------//
// CopyFieldList -- An implementation of UpdatePolicyBase appropriate for
// copying one state field to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CopyFieldList_hh__
#define __Spheral_CopyFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class CopyFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  CopyFieldList(const std::string& masterState, const std::string& copyState);
  virtual ~CopyFieldList();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::string mMasterStateName;
  std::string mCopyStateName;

  CopyFieldList(const CopyFieldList& rhs);
  CopyFieldList& operator=(const CopyFieldList& rhs);
};

}

#endif
