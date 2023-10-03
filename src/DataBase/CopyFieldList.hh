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
#include "CopyState.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class CopyFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  CopyFieldList(const std::string& masterState, const std::string& copyState);
  virtual ~CopyFieldList();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return CopyState<Dimension, ValueType>::prefix(); }

private:
  //--------------------------- Private Interface ---------------------------//
  CopyFieldList(const CopyFieldList& rhs);
  CopyFieldList& operator=(const CopyFieldList& rhs);
};

}

#include "CopyFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class CopyFieldList;
}

#endif
