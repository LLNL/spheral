//---------------------------------Spheral++----------------------------------//
// MassFluxPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//

#include "MassFluxPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
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

  // state
  auto m = state.field(key, 0.0);

  // deriv
  const auto dmdt = derivs.field(this->prefix() + key, 0.0);

// Loop over the internal values of the field.
  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
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
  const auto* rhsPtr = dynamic_cast<const MassFluxPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

