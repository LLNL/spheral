//---------------------------------Spheral++----------------------------------//
// MeshPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the Mesh in the state.
//
// Created by JMO, Sat Feb 12 14:37:57 PST 2011
//----------------------------------------------------------------------------//
#include "MeshPolicy.hh"
#include "Physics/Physics.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using PhysicsSpace::Physics;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor without specifying bounds.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
MeshPolicy(const PhysicsSpace::Physics<Dimension>& package,
           const double voidThreshold):
  UpdatePolicyBase<Dimension>(HydroFieldNames::position + 
                              UpdatePolicyBase<Dimension>::wildcard()),
  mPackage(package),
  mVoidThreshold(voidThreshold),
  mComputeBounds(true),
  mXmin(),
  mXmax() {
}

//------------------------------------------------------------------------------
// Constructor where we specify bounds.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
MeshPolicy(const PhysicsSpace::Physics<Dimension>& package,
           const Vector& xmin,
           const Vector& xmax,
           const double voidThreshold):
  UpdatePolicyBase<Dimension>(HydroFieldNames::position + 
                              UpdatePolicyBase<Dimension>::wildcard()),
  mPackage(package),
  mVoidThreshold(voidThreshold),
  mComputeBounds(false),
  mXmin(xmin),
  mXmax(xmax) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
~MeshPolicy() {
}

//------------------------------------------------------------------------------
// Update the Mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MeshPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  REQUIRE(key == HydroFieldNames::mesh);

  // Find the global bounding box.
  if (mComputeBounds) {
    const FieldSpace::FieldList<Dimension, Vector> positions = state.fields(HydroFieldNames::position, Vector::zero);
    globalBoundingBox<Dimension>(positions, mXmin, mXmax, 
                                 false);     // ghost points
  }

  // This is a special case -- the state knows how to generate the mesh.
  state.generateMesh(mXmin,                     // xmin
                     mXmax,                     // xmax
                     false,                     // generate void
                     false,                     // parallel connectivity
                     mVoidThreshold,
                     mPackage.boundaryBegin(),
                     mPackage.boundaryEnd());
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MeshPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const MeshPolicy<Dimension>* rhsPtr = dynamic_cast<const MeshPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

