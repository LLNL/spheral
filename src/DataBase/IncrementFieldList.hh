//---------------------------------Spheral++----------------------------------//
// IncrementFieldList -- An implementation of FieldListUpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementFieldList_hh__
#define __Spheral_IncrementFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"
#include "IncrementState.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  explicit IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                              const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1, const std::string& depend2,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const FieldList<Dimension, ValueType>& fieldList,
                     const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5, const std::string& depend6
                     const bool wildCardDerivs = false);
  virtual ~IncrementFieldList();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return IncrementState<Dimension, Value>::prefix(); }

private:
  //--------------------------- Private Interface ---------------------------//
  IncrementFieldList(const IncrementFieldList& rhs);
  IncrementFieldList& operator=(const IncrementFieldList& rhs);
};

}

#include "IncrementFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class IncrementFieldList;
}

#endif
