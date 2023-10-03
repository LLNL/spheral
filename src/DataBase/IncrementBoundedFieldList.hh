//---------------------------------Spheral++----------------------------------//
// IncrementBoundedFieldList -- An implementation of FieldListUpdatePolicyBase 
// appropriate for when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementBoundedFieldList_hh__
#define __Spheral_IncrementBoundedFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"
#include "IncrementBoundedState.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class IncrementBoundedFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0, const std::string& depend1,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0, const std::string& depend1, const std::string& depend2,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                            const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                            const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                            const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  virtual ~IncrementBoundedFieldList();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return IncrementBoundedState<Dimension, ValueType>::prefix(); }

private:
  //--------------------------- Private Interface ---------------------------//
  IncrementBoundedFieldList(const IncrementBoundedFieldList& rhs);
  IncrementBoundedFieldList& operator=(const IncrementBoundedFieldList& rhs);
};

}

#include "IncrementBoundedFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename DataType, typename BoundValueType> class IncrementBoundedFieldList;
}

#endif
