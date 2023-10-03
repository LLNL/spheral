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
  IncrementState(const bool wildCardDerivs = false);
  explicit IncrementState(const std::string& depend0,
                          const bool wildCardDerivs = false);
  IncrementState(const std::string& depend0, const std::string& depend1,
                 const bool wildCardDerivs = false);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2,
                 const bool wildCardDerivs = false);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                 const bool wildCardDerivs = false);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                 const bool wildCardDerivs = false);
  IncrementState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                 const bool wildCardDerivs = false);
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

  // Flip whether we try to find multiple registered increment fields.
  bool wildCardDerivs() const;
  void wildCardDerivs(const bool val);

private:
  //--------------------------- Private Interface ---------------------------//
  IncrementState(const IncrementState& rhs);
  IncrementState& operator=(const IncrementState& rhs);

  // Flag for looking for multiple increment derivatives.
  bool mWildCardDerivs;
};

}

#include "IncrementStateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class IncrementState;
}

#endif
