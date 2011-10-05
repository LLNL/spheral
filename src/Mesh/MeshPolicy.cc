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
#include "Utilities/boundingVolumes.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using PhysicsSpace::Physics;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
MeshPolicy(const PhysicsSpace::Physics<Dimension>& package,
           const double voidThreshold):
  UpdatePolicyBase<Dimension>(HydroFieldNames::position + 
                              UpdatePolicyBase<Dimension>::wildcard()),
  mPackage(package),
  mVoidThreshold(voidThreshold) {
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
  Vector xmin, xmax;
  const FieldSpace::FieldList<Dimension, Vector> positions = state.fields(HydroFieldNames::position, Vector::zero);
  boundingBox<Dimension>(positions, xmin, xmax, 
                         false,      // ghost points
                         true);      // quantize results

  // Puff things up a bit.
  const Vector delta = xmax - xmin;
  xmin -= 0.01*delta;
  xmax += 0.01*delta;

  // This is a special case -- the state knows how to generate the mesh.
  state.generateMesh(xmin,                     // xmin
                     xmax,                     // xmax
                     false,                    // generate void
                     false,                    // parallel connectivity
                     mVoidThreshold,
                     mPackage.boundaryBegin(),
                     mPackage.boundaryEnd());

//   // Blago!
//   double xmax = -1e10;
//   const MeshSpace::Mesh<Dimension>& mesh = state.mesh();
//   for (unsigned i = 0; i != mesh.numNodes(); ++i) xmax = std::max(xmax, mesh.node(i).position().x());
//   cerr << "Max node position:  " << xmax << endl;
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

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class MeshPolicy<Dim<1> >;
  template class MeshPolicy<Dim<2> >;
  template class MeshPolicy<Dim<3> >;
}

