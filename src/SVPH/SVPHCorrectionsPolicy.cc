//---------------------------------Spheral++----------------------------------//
// SVPHCorrectionsPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the SVPHCorrections in the state.
//
// Created by JMO, Fri Aug  9 15:24:04 PDT 2013
//----------------------------------------------------------------------------//
#include "SVPHCorrectionsPolicy.hh"
#include "computeSVPHCorrections.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using DataBaseSpace::DataBase;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHCorrectionsPolicy<Dimension>::
SVPHCorrectionsPolicy(const DataBase<Dimension>& dataBase,
                      const TableKernel<Dimension>& kernel):
  UpdatePolicyBase<Dimension>(HydroFieldNames::position + UpdatePolicyBase<Dimension>::wildcard(),
                              HydroFieldNames::H,
                              HydroFieldNames::volume),
  mDataBase(dataBase),
  mKernel(kernel) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHCorrectionsPolicy<Dimension>::
~SVPHCorrectionsPolicy() {
}

//------------------------------------------------------------------------------
// Update the SVPHCorrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHCorrectionsPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  const KeyType Akey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::A_CSPH, nodeListKey);
  const KeyType Bkey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::B_CSPH, nodeListKey);
  const KeyType gradBkey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::gradB_CSPH, nodeListKey);

  Field<Dimension, Scalar>& A = state.field(Akey, 0.0);
  Field<Dimension, Vector>& B = state.field(Bkey, Vector::zero);
  Field<Dimension, Tensor>& gradB = state.field(gradBkey, Tensor::zero);

  SVPHSpace::computeSVPHCorrections<Dimension>(mDataBase.connectivityMap(),
                                               mKernel,
                                               volume, position, H,
                                               A, B, gradB);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SVPHCorrectionsPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a SVPHCorrections object.
  const SVPHCorrectionsPolicy<Dimension>* rhsPtr = dynamic_cast<const SVPHCorrectionsPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

