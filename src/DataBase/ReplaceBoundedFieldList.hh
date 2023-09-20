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

#include <float.h>
#include "FieldListUpdatePolicyBase.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class ReplaceBoundedFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  ReplaceBoundedFieldList(const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0, const std::string& depend1,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                          const BoundValueType minValue = BoundValueType(-DBL_MAX),
                          const BoundValueType maxValue = BoundValueType(DBL_MAX));
  virtual ~ReplaceBoundedFieldList();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // An alternate method to be called when you want to specify that the "Replace" information
  // in the derivatives is invalid, and instead the value should be treated as a time advancement
  // algorithm instead.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt);

  // Access the min and max's.
  BoundValueType minValue() const;
  BoundValueType maxValue() const;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "new "; }


private:
  //--------------------------- Private Interface ---------------------------//
  BoundValueType mMinValue;
  BoundValueType mMaxValue;

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
