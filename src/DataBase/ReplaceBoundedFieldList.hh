//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceBoundedFieldList_hh__
#define __Spheral_ReplaceBoundedFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"
#include "ReplaceBoundedState.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class ReplaceBoundedFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  explicit ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                                   const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                                   const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0, const std::string& depend1,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0, const std::string& depend1, const std::string& depend2,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  virtual ~ReplaceBoundedFieldList();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return ReplaceBoundedState<Dimension, ValueType>::prefix(); }

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceBoundedFieldList(const ReplaceBoundedFieldList& rhs);
  ReplaceBoundedFieldList& operator=(const ReplaceBoundedFieldList& rhs);
};

}

#include "ReplaceBoundedFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename DataType, typename BoundValueType> class ReplaceBoundedFieldList;
}

#endif
