//---------------------------------Spheral++----------------------------------//
// VoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density according to the specific volume from the
// Voronoi tesselation.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "VoronoiMassDensityPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiMassDensityPolicy<Dimension>::
VoronoiMassDensityPolicy(const double rhoMin, const double rhoMax):
  ReplaceState<Dimension, typename Dimension::Scalar>(HydroFieldNames::mass,
                                                      HydroFieldNames::volume),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
  REQUIRE(rhoMin <= rhoMax);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiMassDensityPolicy<Dimension>::
~VoronoiMassDensityPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity);
  Field<Dimension, Scalar>& massDensity = state.field(key, 0.0);

  // Get the mass and volume from the state.
  const KeyType massKey = State<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const KeyType volKey = State<Dimension>::buildFieldKey(HydroFieldNames::volume, nodeListKey);
  CHECK(state.registered(massKey));
  CHECK(state.registered(volKey));
  const Field<Dimension, Scalar>& mass = state.field(massKey, 0.0);
  const Field<Dimension, Scalar>& volume = state.field(volKey, 0.0);

  // Set the mass density.
  const NodeList<Dimension>& nodeList = mass.nodeList();
  for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) {
    massDensity(i) = max(mRhoMin, min(mRhoMax, mass(i) * safeInv(volume(i))));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VoronoiMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const VoronoiMassDensityPolicy<Dimension>* rhsPtr = dynamic_cast<const VoronoiMassDensityPolicy<Dimension>*>(&rhs);
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
  template class VoronoiMassDensityPolicy<Dim<1> >;
  template class VoronoiMassDensityPolicy<Dim<2> >;
  template class VoronoiMassDensityPolicy<Dim<3> >;
}

