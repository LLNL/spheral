//---------------------------------Spheral++----------------------------------//
// IncrementBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementBoundedState_hh__
#define __Spheral_IncrementBoundedState_hh__

#include <float.h>
#include "FieldUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class IncrementBoundedState: public FieldUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  IncrementBoundedState(const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0, const std::string& depend1,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  IncrementBoundedState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                        const BoundValueType minValue = BoundValueType(-DBL_MAX),
                        const BoundValueType maxValue = BoundValueType(DBL_MAX));
  virtual ~IncrementBoundedState();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
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

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  BoundValueType mMinValue;
  BoundValueType mMaxValue;

  IncrementBoundedState(const IncrementBoundedState& rhs);
  IncrementBoundedState& operator=(const IncrementBoundedState& rhs);
};

}

#ifndef __GCCXML__
#include "IncrementBoundedStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename DataType, typename BoundValueType> class IncrementBoundedState;
}

#endif
