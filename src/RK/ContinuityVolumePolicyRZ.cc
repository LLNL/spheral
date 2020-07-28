//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Specialized for RZ geometry.
//
// Created by JMO, Tue Oct 29 16:14:03 PDT 2019
//----------------------------------------------------------------------------//

#include "ContinuityVolumePolicyRZ.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


namespace {
//------------------------------------------------------------------------------
// Compute the dimension dependent volume of the H tensor.
//------------------------------------------------------------------------------
double Hvolume(const Dim<2>::SymTensor& H) {
  return M_PI/H.Determinant();
}
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
ContinuityVolumePolicyRZ::
ContinuityVolumePolicyRZ():
  IncrementFieldList<Dim<2>, Dim<2>::Scalar>(HydroFieldNames::mass,
                                             HydroFieldNames::massDensity) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
ContinuityVolumePolicyRZ::
~ContinuityVolumePolicyRZ() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
ContinuityVolumePolicyRZ::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  SPHERAL_UNUSED(t);
  SPHERAL_UNUSED(dt);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> volume = state.fields(fieldKey, Scalar());
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto DrhoDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);

  // Loop over the internal values of the field.
  const auto numNodeLists = volume.size();
  for (auto k = 0u; k != numNodeLists; ++k) {
    const auto n = volume[k]->numInternalElements();
    for (auto i = 0u; i != n; ++i) {
      const auto circi = 2.0*M_PI*abs(pos(k,i).y());
      CHECK(circi > 0.0);
      const auto volMin = 0.5*mass(k,i)*safeInvVar(circi*rho(k,i));
      const auto volMax = Hvolume(H(k,i));
      const auto dVdt = -mass(k,i)*safeInvVar(circi*rho(k,i)*rho(k,i))*DrhoDt(k,i);
      volume(k,i) = std::max(volMin, std::min(volMax, volume(k,i) + multiplier*dVdt));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
ContinuityVolumePolicyRZ::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ContinuityVolumePolicyRZ* rhsPtr = dynamic_cast<const ContinuityVolumePolicyRZ*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

