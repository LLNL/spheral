//---------------------------------Spheral++----------------------------------//
// DeviatoricStressPolity -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#include "DeviatoricStressPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Kernel/TableKernel.hh"
#include "Material/EquationOfState.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
DeviatoricStressPolicy<Dimension>::
DeviatoricStressPolicy():
  FieldUpdatePolicy<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DeviatoricStressPolicy<Dimension>::
~DeviatoricStressPolicy() {
}

//------------------------------------------------------------------------------
// Update the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DeviatoricStressPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::deviatoricStress);

  // Get the state we're advancing.
  auto&       S = state.field(key, SymTensor::zero);
  const auto& DSDt = derivs.field(IncrementState<Dimension, SymTensor>::prefix() + key, SymTensor::zero);

  // We only want to enforce zeroing the trace in Cartesian coordinates.   In RZ or R
  // we assume the missing components on the diagonal sum to -Trace(S).
  const auto zeroTrace = GeometryRegistrar::coords() == CoordinateType::Cartesian;

  // Iterate over the internal nodes.
  const auto n = S.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    auto S0 = S(i) + multiplier*(DSDt(i));               // Elastic prediction for the new deviatoric stress
    if (zeroTrace) {
      S0 -= SymTensor::one * S0.Trace()/Dimension::nDim; // Ensure the deviatoric stress is traceless (all but RZ and spherical)
      CHECK(fuzzyEqual(S0.Trace(), 0.0));
    }

    // Purely elastic flow.  The plastic yielding is accounted for when we update the plastic strain.
    S(i) = S0;
  }

//     // Finally apply the pressure limits to the allowed deviatoric stress.
//     S(i) = max(Pmin, min(Pmax, S(i)));
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DeviatoricStressPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is a DeviatoricStress operator, and has
  // the same cutoff values.
  const auto rhsPtr = dynamic_cast<const DeviatoricStressPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

