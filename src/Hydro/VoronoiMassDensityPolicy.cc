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
#include "Field/FieldList.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiMassDensityPolicy<Dimension>::
VoronoiMassDensityPolicy(const double rhoMin, const double rhoMax):
  ReplaceFieldList<Dimension, typename Dimension::Scalar>(HydroFieldNames::mass,
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
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> massDensity = state.fields(fieldKey, Scalar());
  const unsigned numFields = massDensity.numFields();

  // Get the mass and volume from the state.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);

  // Set the mass density.
  for (unsigned i = 0; i != numFields; ++i) {
    for (unsigned j = 0; j != massDensity[i]->numInternalElements(); ++j) {
      massDensity(i,j) = max(mRhoMin, min(mRhoMax, mass(i,j) * safeInv(volume(i,j))));
    }
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

