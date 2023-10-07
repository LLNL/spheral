//---------------------------------Spheral++----------------------------------//
// CopyState -- An implementation of UpdatePolicyBase appropriate for
// copying one state field to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Tue Oct 5 11:08:48 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_CopyState_hh__
#define __Spheral_CopyState_hh__

#include "UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType>
class CopyState: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  CopyState(const std::string& masterState, const std::string& copyState);
  virtual ~CopyState();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::string mMasterStateName;
  std::string mCopyStateName;

  CopyState(const CopyState& rhs);
  CopyState& operator=(const CopyState& rhs);
};

}

#include "CopyStateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class CopyState;
}

#endif
