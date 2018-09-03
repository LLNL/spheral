//---------------------------------Spheral++----------------------------------//
// RZDeviatoricStressPolicy.
//
// Specialized version of the IncrementFieldList policy, with some criteria for
// zeroing out the deviatoric stress in special cases.
// This one is also specialized for the RZ strength case.
//
// Created by JMO, Mon May  9 14:09:12 PDT 2016
//----------------------------------------------------------------------------//
#include "RZDeviatoricStressPolicy.hh"
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

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
RZDeviatoricStressPolicy::
RZDeviatoricStressPolicy():
  IncrementFieldList<Dim<2>, Dim<2>::SymTensor>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
RZDeviatoricStressPolicy::
~RZDeviatoricStressPolicy() {
}

//------------------------------------------------------------------------------
// Update the FieldList.
//------------------------------------------------------------------------------
void
RZDeviatoricStressPolicy::
update(const KeyType& key,
       State<Dim<2> >& state,
       StateDerivatives<Dim<2> >& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the state we're advancing.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  const FieldList<Dimension, SymTensor> DSDt = derivs.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + 
                                                             SolidFieldNames::deviatoricStress, SymTensor::zero);
  const FieldList<Dimension, Scalar> DSTTDt = derivs.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + 
                                                            SolidFieldNames::deviatoricStressTT, 0.0);

  // Iterate over the internal nodes.
  const unsigned numFields = S.numFields();
  for (unsigned k = 0; k != numFields; ++k) {
    const unsigned n = S[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const SymTensor Sold = S(k,i);                          // Starting deviatoric stress.
      const SymTensor S0 = Sold + multiplier*(DSDt(k,i));     // Elastic prediction for the new deviatoric stress.

      const Scalar STTold = STT(k,i);
      const Scalar STT0 = STTold + multiplier*(DSTTDt(k,i));

      // Purely elastic flow.  The plastic yielding is accounted for when we update the plastic strain.
      S(k,i) = S0;
      STT(k,i) = STT0;
    }
  }

//     // Finally apply the pressure limits to the allowed deviatoric stress.
//     S(i) = max(Pmin, min(Pmax, S(i)));
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
RZDeviatoricStressPolicy::
operator==(const UpdatePolicyBase<Dim<2> >& rhs) const {

  // We're only equal if the other guy is a DeviatoricStress operator, and has
  // the same cutoff values.
  const RZDeviatoricStressPolicy* rhsPtr = dynamic_cast<const RZDeviatoricStressPolicy*>(&rhs);
  return (rhsPtr != 0);
}

}

