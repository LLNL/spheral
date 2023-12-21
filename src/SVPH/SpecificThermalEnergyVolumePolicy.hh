//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyVolumePolicy
//
// Created by JMO, Sun Aug 18 19:21:27 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SpecificThermalEnergyVolumePolicy_hh__
#define __Spheral_SpecificThermalEnergyVolumePolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

template<typename Dimension>
class SpecificThermalEnergyVolumePolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  SpecificThermalEnergyVolumePolicy();
  virtual ~SpecificThermalEnergyVolumePolicy();
  
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
  SpecificThermalEnergyVolumePolicy(const SpecificThermalEnergyVolumePolicy& rhs);
  SpecificThermalEnergyVolumePolicy& operator=(const SpecificThermalEnergyVolumePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SpecificThermalEnergyVolumePolicy;
}

#endif
