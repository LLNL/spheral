//---------------------------------Spheral++----------------------------------//
// MassFluxPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//

#include "MassFluxPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------ 
template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(std::initializer_list<std::string> depends):
  IncrementState<Dimension, Scalar>(depends) {
}
  
//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MassFluxPolicy<Dimension>::
~MassFluxPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MassFluxPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::mass and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // state
  auto m = state.field(fieldKey, 0.0);

  // deriv
  const auto dmdt = derivs.field(IncrementState<Dimension, Scalar>::prefix() + fieldKey, 0.0);

// Loop over the internal values of the field.
  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (auto j = 0u; j < n; ++j) {
    m(i) += std::max(multiplier*dmdt(i),-m(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MassFluxPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const MassFluxPolicy<Dimension>* rhsPtr = dynamic_cast<const MassFluxPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

