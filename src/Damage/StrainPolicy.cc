//---------------------------------Spheral++----------------------------------//
// StrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the strain.
//
// Created by JMO, Tue Oct 12 16:43:23 PDT 2004
//----------------------------------------------------------------------------//
#include "StrainPolicy.hh"
#include "NodeList/NodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Kernel/HatKernel.hh"

#include <vector>
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
StrainPolicy<Dimension>::
StrainPolicy():
  UpdatePolicyBase<Dimension>({HydroFieldNames::position,
                               HydroFieldNames::H,
                               SolidFieldNames::YoungsModulus,
                               HydroFieldNames::pressure,
                               SolidFieldNames::deviatoricStress}) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrainPolicy<Dimension>::
~StrainPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::strain);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  const double tiny = 1.0e-30;

  // Get the state fields.
  const KeyType EKey = State<Dimension>::buildFieldKey(SolidFieldNames::YoungsModulus, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType stressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  CHECK(state.registered(EKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(stressKey));

  const Field<Dimension, Scalar>& E = state.field(EKey, 0.0);
  const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);
  const Field<Dimension, SymTensor>& S = state.field(stressKey, SymTensor::zero);

  // Iterate over the internal nodes.
  for (auto i = 0u; i != stateField.numInternalElements(); ++i) {

    // Compute the maximum tensile stress.
//     Scalar Pi = P(i);
//     if (Pi < 0.0) Pi *= 1.0 - D(i);
    const SymTensor sigmai = S(i) - P(i)*SymTensor::one;
    const Scalar sigmati = max(0.0, sigmai.eigenValues().maxElement());
    CHECK(sigmati >= 0.0);

    // Compute the initial strain for this node.
    CHECK(E(i) >= 0.0);
    stateField(i) = sigmati/(E(i) + tiny*max(1.0, sigmati));
    CHECK(stateField(i) >= 0.0);

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const StrainPolicy<Dimension>* rhsPtr = dynamic_cast<const StrainPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

