//---------------------------------Spheral++----------------------------------//
// HVolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the local hull constructions.
//
// Created by JMO, Wed Aug 13 10:52:16 PDT 2014
//----------------------------------------------------------------------------//

#include "HVolumePolicy.hh"
#include "computeHVolumes.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::FieldList;
using MeshSpace::Mesh;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HVolumePolicy<Dimension>::
HVolumePolicy(const Scalar kernelExtent):
  mKernelExtent(kernelExtent),
  ReplaceFieldList<Dimension, typename Dimension::Scalar>(HydroFieldNames::H) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HVolumePolicy<Dimension>::
~HVolumePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HVolumePolicy<Dimension>::
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

  // Get the H field from the state, and do the deed.
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  computeHVolumes(mKernelExtent, H, volume);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const HVolumePolicy<Dimension>* rhsPtr = dynamic_cast<const HVolumePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

