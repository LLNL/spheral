//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Created by JMO, Tue Sep 20 14:53:32 PDT 2016
//----------------------------------------------------------------------------//

#include "ContinuityVolumePolicy.hh"
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
double Hvolume(const Dim<1>::SymTensor& H) {
  return 2.0/H.xx();
}

double Hvolume(const Dim<2>::SymTensor& H) {
  return M_PI/H.Determinant();
}

double Hvolume(const Dim<3>::SymTensor& H) {
  return 4.0*M_PI/(3.0*H.Determinant());
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ContinuityVolumePolicy<Dimension>::
ContinuityVolumePolicy():
  IncrementFieldList<Dimension, typename Dimension::Scalar>(HydroFieldNames::mass,
                                                            HydroFieldNames::massDensity) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ContinuityVolumePolicy<Dimension>::
~ContinuityVolumePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ContinuityVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> volume = state.fields(fieldKey, Scalar());
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> DrhoDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);

  // Loop over the internal values of the field.
  const unsigned numNodeLists = volume.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = volume[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar volMin = 0.5*mass(k,i)*safeInvVar(rho(k,i));
      const Scalar volMax = Hvolume(H(k,i));
      const Scalar dVdt = -mass(k,i)*safeInvVar(rho(k,i)*rho(k,i))*DrhoDt(k,i);
      volume(k,i) = std::max(volMin, std::min(volMax, volume(k,i) + multiplier*dVdt));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ContinuityVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ContinuityVolumePolicy<Dimension>* rhsPtr = dynamic_cast<const ContinuityVolumePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

