//---------------------------------Spheral++----------------------------------//
// SpecificFromTotalThermalEnergyPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the specific thermal energy from the total
// energy.
// 
// Created by JMO, Thu Jan  7 23:08:39 PST 2016
//----------------------------------------------------------------------------//

#include <vector>

#include "SpecificFromTotalThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificFromTotalThermalEnergyPolicy<Dimension>::
SpecificFromTotalThermalEnergyPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificFromTotalThermalEnergyPolicy<Dimension>::
~SpecificFromTotalThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificFromTotalThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // Get the state fields.
  const auto mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto DEDt = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);

  // Do it.
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {
    const size_t n = eps[nodeListi]->numInternalElements();
    for (size_t i = 0; i != n; ++i) {
      auto& epsi = eps(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi0 = velocity(nodeListi, i);
      const auto& ai0 = DvDt(nodeListi, i);
      const auto  DEDti = DEDt(nodeListi, i);
      const auto  E0i = mi*(0.5*vi0.magnitude2() + epsi);
      const auto  E1i = E0i + multiplier*DEDti;
      epsi = E1i/mi - 0.5*(vi0 + multiplier*ai0).magnitude2();
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpecificFromTotalThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const SpecificFromTotalThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const SpecificFromTotalThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

