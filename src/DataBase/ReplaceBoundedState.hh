//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Tue Aug 31 14:03:45 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceBoundedState_hh__
#define __Spheral_ReplaceBoundedState_hh__

#include <float.h>
#include "FieldUpdatePolicyBase.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class ReplaceBoundedState: public FieldUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  ReplaceBoundedState(const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0, const std::string& depend1,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  ReplaceBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                      const BoundValueType minValue = BoundValueType(-DBL_MAX),
                      const BoundValueType maxValue = BoundValueType(DBL_MAX));
  virtual ~ReplaceBoundedState();
  
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

  ReplaceBoundedState(const ReplaceBoundedState& rhs);
  ReplaceBoundedState& operator=(const ReplaceBoundedState& rhs);
};

}

#ifndef __GCCXML__
#include "ReplaceBoundedStateInline.hh"
#endif

#endif
