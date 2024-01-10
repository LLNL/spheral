//---------------------------------Spheral++----------------------------------//
// MFVIncrementSpecificThermalEnergyPolicy -- policy to update the velocity from the 
//                                  momentum time derivative
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral_MFVIncrementSpecificThermalEnergyPolicy_hh__
#define __Spheral_MFVIncrementSpecificThermalEnergyPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

template<typename Dimension>
class MFVIncrementSpecificThermalEnergyPolicy: public FieldUpdatePolicy<Dimension> {
public:

  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  MFVIncrementSpecificThermalEnergyPolicy(std::initializer_list<std::string> depends={});
  ~MFVIncrementSpecificThermalEnergyPolicy();
  
  // Overload the methods describing how to update FieldLists.
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

  const std::string mStateKey;
  const std::string mDerivativeKey;

  //--------------------------- Private Interface ---------------------------//
  MFVIncrementSpecificThermalEnergyPolicy(const MFVIncrementSpecificThermalEnergyPolicy& rhs);
  MFVIncrementSpecificThermalEnergyPolicy& operator=(const MFVIncrementSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MFVIncrementSpecificThermalEnergyPolicy;
}

#endif
