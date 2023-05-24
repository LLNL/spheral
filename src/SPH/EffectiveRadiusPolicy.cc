//---------------------------------Spheral++----------------------------------//
// EffectiveRadiusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the grad h correction terms.
//
// Created by JMO, Tue Oct 30 15:42:41 PDT 2007
//----------------------------------------------------------------------------//
#include "SPH/EffectiveRadiusPolicy.hh"
#include "SPH/SPHHydroBaseRZ.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
EffectiveRadiusPolicy::
EffectiveRadiusPolicy():
  FieldListUpdatePolicyBase<Dimension, Scalar>(HydroFieldNames::position,
                                               HydroFieldNames::H) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
EffectiveRadiusPolicy::
~EffectiveRadiusPolicy() {
}

//------------------------------------------------------------------------------
// Update the FieldList.
//------------------------------------------------------------------------------
void
EffectiveRadiusPolicy::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  using Vector = Dim<2>::Vector;
  using SymTensor = Dim<2>::SymTensor;

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  REQUIRE(fieldKey == HydroFieldNames::reff);

  // Grab the state we need.
  auto reff = state.fields(fieldKey, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);

  // Walk all the points.
  const auto numFields = reff.numFields();
  CHECK(position.numFields() == numFields and
        H.numFields() == numFields);
  for (auto k = 0u; k < numFields; ++k) {
    const auto& nodeList = reff[k]->nodeList();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();
    const auto  n = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& posi = position(k, i);
      const auto& Hi = H(k, i);
      const auto  zetai = abs((Hi*posi).y());
      const auto  hri = abs(posi.y())*safeInv(zetai);
      reff(k, i) = SPHHydroBaseRZ::reff(posi.y(), hri, nPerh);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
EffectiveRadiusPolicy::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const EffectiveRadiusPolicy* rhsPtr = dynamic_cast<const EffectiveRadiusPolicy*>(&rhs);
  if (rhsPtr == nullptr) {
    return false;
  } else {
    return true;
  }
}

}
