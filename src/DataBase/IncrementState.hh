//---------------------------------Spheral++----------------------------------//
// IncrementState -- An implementation of FieldUpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementState_hh__
#define __Spheral_IncrementState_hh__

#include "FieldUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementState: public FieldUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  IncrementState();
  explicit IncrementState(const std::string& depend0);
  IncrementState(const std::string& depend0, const std::string& depend1);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~IncrementState();
  
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
  IncrementState(const IncrementState& rhs);
  IncrementState& operator=(const IncrementState& rhs);
};

}

#include "IncrementStateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class IncrementState;
}

#endif
