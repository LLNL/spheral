//---------------------------------Spheral++----------------------------------//
// DamagedPressurePolicy -- Override the default sound speed policy in the 
// presence of damage.
//
// Created by JMO, Thu Jan 28 14:23:46 PST 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamagedPressurePolicy_hh__
#define __Spheral_DamagedPressurePolicy_hh__

#include "Hydro/PressurePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class DamagedPressurePolicy: public PressurePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  DamagedPressurePolicy();
  virtual ~DamagedPressurePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  DamagedPressurePolicy(const DamagedPressurePolicy& rhs);
  DamagedPressurePolicy& operator=(const DamagedPressurePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DamagedPressurePolicy;
}

#endif
