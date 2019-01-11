//---------------------------------Spheral++----------------------------------//
// MashCorrectionPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the MASH node based corrections.
//
// Created by JMO, Mon Sep 26 21:23:07 PDT 2005
//----------------------------------------------------------------------------//

#include "MashCorrectionPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/MashNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MashCorrectionPolicy<Dimension>::
MashCorrectionPolicy():
  UpdatePolicyBase<Dimension, FieldType>(HydroFieldNames::position,
                                         HydroFieldNames::H,
                                         HydroFieldNames::weight) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MashCorrectionPolicy<Dimension>::
~MashCorrectionPolicy() {
}

//------------------------------------------------------------------------------
// Update the MASH normalization and gradient corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MashCorrectionPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  REQUIRE(key.second == HydroFieldNames::mashCorrection);

  typedef typename Dimension::Scalar Scalar;

  // Get the normalization and correction fields.
  const KeyType normalizationKey(key.first, HydroFieldNames::mashNormalization);
  const KeyType correctionKey(key.first, HydroFieldNames::mashCorrection);
  const KeyType QnormalizationKey(key.first, HydroFieldNames::mashQnormalization);
  const KeyType QcorrectionKey(key.first, HydroFieldNames::mashQcorrection);
  CHECK(state.scalarFieldRegistered(normalizationKey));
  CHECK(state.tensorFieldRegistered(correctionKey));
  CHECK(state.scalarFieldRegistered(QnormalizationKey));
  CHECK(state.tensorFieldRegistered(QcorrectionKey));
  Field<Dimension, Scalar>& normalization = state.scalarField(normalizationKey);
  Field<Dimension, Tensor>& correction = state.tensorField(correctionKey);
  Field<Dimension, Scalar>& Qnormalization = state.scalarField(QnormalizationKey);
  Field<Dimension, Tensor>& Qcorrection = state.tensorField(QcorrectionKey);

  // We use the virtual method of the MashNodeList to do the deed.
  const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(key.first);
  CHECK(fluidNodeListPtr != 0);
  fluidNodeListPtr->computeCorrections(normalization, correction, 
                                       Qnormalization, Qcorrection,
                                       state);

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MashCorrectionPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension, FieldType>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const MashCorrectionPolicy<Dimension>* rhsPtr = dynamic_cast<const MashCorrectionPolicy<Dimension>*>(&rhs);
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
  template class MashCorrectionPolicy<Dim<1> >;
  template class MashCorrectionPolicy<Dim<2> >;
  template class MashCorrectionPolicy<Dim<3> >;
}
