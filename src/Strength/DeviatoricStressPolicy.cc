//---------------------------------Spheral++----------------------------------//
// DeviatoricStressPolity -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#include "DeviatoricStressPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Kernel/TableKernel.hh"
#include "Material/EquationOfState.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
DeviatoricStressPolicy<Dimension>::
DeviatoricStressPolicy():
  IncrementFieldList<Dimension, typename Dimension::SymTensor>() {
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

  // Get the state we're advancing.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const FieldList<Dimension, SymTensor> DSDt = derivs.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + 
                                                             SolidFieldNames::deviatoricStress, SymTensor::zero);

  // Iterate over the internal nodes.
  const unsigned numFields = S.numFields();
  for (unsigned k = 0; k != numFields; ++k) {
    const unsigned n = S[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const SymTensor Sold = S(k,i);                         // Starting deviatoric stress.
      const SymTensor S0 = Sold + multiplier*(DSDt(k,i));    // Elastic prediction for the new deviatoric stress.

      // Purely elastic flow.  The plastic yielding is accounted for when we update the plastic strain.
      S(k,i) = S0;
    }
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
  const DeviatoricStressPolicy<Dimension>* rhsPtr = dynamic_cast<const DeviatoricStressPolicy<Dimension>*>(&rhs);
  return (rhsPtr != 0);
}

}

