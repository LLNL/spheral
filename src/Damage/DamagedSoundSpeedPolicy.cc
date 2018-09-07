//---------------------------------Spheral++----------------------------------//
// DamagedSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of damage.
//----------------------------------------------------------------------------//
#include "DamagedSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamagedSoundSpeedPolicy<Dimension>::
DamagedSoundSpeedPolicy():
  SoundSpeedPolicy<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamagedSoundSpeedPolicy<Dimension>::
~DamagedSoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamagedSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);

  // Get the sound speed.
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // Have the base class set the initial pressure.
  SoundSpeedPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);

  // Is there a scalar damage field for this NodeList?
  {
    const KeyType DKey = State<Dimension>::buildFieldKey(SolidFieldNames::scalarDamage, nodeListKey);
    if (state.registered(DKey)) {
      const Field<Dimension, Scalar>& D = state.field(DKey, 0.0);
      for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
        CHECK(D(i) >= 0.0 && D(i) <= 1.0);
        const Scalar fDi = min(0.0, max(1.0, 1.0 - D(i)));
        stateField(i) *= fDi*fDi;
      }
    }
  }

  // Is there a tensor damage field for this NodeList?
  {
    typedef typename Dimension::SymTensor SymTensor;
    const KeyType DKey = State<Dimension>::buildFieldKey(SolidFieldNames::effectiveTensorDamage, nodeListKey);
    if (state.registered(DKey)) {
      const Field<Dimension, SymTensor>& D = state.field(DKey, SymTensor::zero);
      for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
        const Scalar fDi = max(0.0, min(1.0, 1.0 - D(i).eigenValues().maxElement()));
        CHECK(fDi >= 0.0 && fDi <= 1.0);
        stateField(i) *= fDi*fDi;
      }
    }
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DamagedSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const DamagedSoundSpeedPolicy<Dimension>* rhsPtr = dynamic_cast<const DamagedSoundSpeedPolicy<Dimension>*>(&rhs);
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
  template class DamagedSoundSpeedPolicy<Dim<1> >;
  template class DamagedSoundSpeedPolicy<Dim<2> >;
  template class DamagedSoundSpeedPolicy<Dim<3> >;
}

