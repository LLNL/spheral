//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyVolumePolicy
//
// Created by JMO, Sun Aug 18 19:21:27 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SpecificThermalEnergyVolumePolicy_hh__
#define __Spheral_SpecificThermalEnergyVolumePolicy_hh__

#include "DataBase/FieldUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension>
class SpecificThermalEnergyVolumePolicy: 
    public FieldUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  SpecificThermalEnergyVolumePolicy();
  virtual ~SpecificThermalEnergyVolumePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  SpecificThermalEnergyVolumePolicy(const SpecificThermalEnergyVolumePolicy& rhs);
  SpecificThermalEnergyVolumePolicy& operator=(const SpecificThermalEnergyVolumePolicy& rhs);
};

}

#endif
