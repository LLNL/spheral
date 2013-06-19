//---------------------------------Spheral++----------------------------------//
// VolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the current Voronoi tesselation.
//
// Created by JMO, Mon Aug  1 15:17:36 PDT 2011
//----------------------------------------------------------------------------//

#include "VolumePolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;
using MeshSpace::Mesh;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VolumePolicy<Dimension>::
VolumePolicy():
  ReplaceState<Dimension, typename Dimension::Scalar>(HydroFieldNames::mesh) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VolumePolicy<Dimension>::
~VolumePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume);
  Field<Dimension, Scalar>& volume = state.field(key, 0.0);

  // Get the mesh from the state.
  const Mesh<Dimension>& mesh = state.mesh();

  // Just read the cell volumes from the mesh as our new value.
  const NodeList<Dimension>& nodeList = volume.nodeList();
  const unsigned offset = mesh.offset(nodeList);
  for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) {
    volume(i) = mesh.zone(offset + i).volume();
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const VolumePolicy<Dimension>* rhsPtr = dynamic_cast<const VolumePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

