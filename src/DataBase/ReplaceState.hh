//---------------------------------Spheral++----------------------------------//
// ReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceState_hh__
#define __Spheral_ReplaceState_hh__

#include "FieldUpdatePolicy.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class ReplaceState: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  ReplaceState();
  explicit ReplaceState(const std::string& depend0);
  ReplaceState(const std::string& depend0, const std::string& depend1);
  ReplaceState(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceState(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceState();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // An alternate method to be called when you want to specify that the "Replace" information
  // in the derivatives is invalid, and instead the value should be treated as a time advancement
  // algorithm instead.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "new "; }

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceState(const ReplaceState& rhs);
  ReplaceState& operator=(const ReplaceState& rhs);
};

}

#include "ReplaceStateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class ReplaceState;
}

#endif
