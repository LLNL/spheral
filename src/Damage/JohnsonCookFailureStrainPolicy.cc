//---------------------------------Spheral++----------------------------------//
// JohnsonCookFailureStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the JC failure strain.
//
// Created by JMO, Thu Jul 12 13:40:45 PDT 2018
//----------------------------------------------------------------------------//

#include <vector>

#include "JohnsonCookFailureStrainPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookFailureStrainPolicy<Dimension>::
JohnsonCookFailureStrainPolicy(const FieldSpace::Field<Dimension, Scalar>& D1,
                               const FieldSpace::Field<Dimension, Scalar>& D2,
                               const double D3,
                               const double D4,
                               const double D5,
                               const double epsilondot0,
                               const double sigmamax,
                               const double efailmin,
                               const double Tcrit):
  UpdatePolicyBase<Dimension>(HydroFieldNames::pressure,
                              HydroFieldNames::specificThermalEnergy,
                              SolidFieldNames::deviatoricStress,
                              SolidFieldNames::plasticStrain,
                              SolidFieldNames::meltSpecificEnergy),
  mD1(D1),
  mD2(D2),
  mD3(D3),
  mD4(D4),
  mD5(D5),
  mepsilondot0(epsilondot0),
  msigmamax(sigmamax),
  mefailmin(efailmin),
  mTcrit(Tcrit) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookFailureStrainPolicy<Dimension>::
~JohnsonCookFailureStrainPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookFailureStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::flaws);
  auto& efail = state.field(key, 0.0);

  // Get the state fields.
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType epsKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  const KeyType stressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  const KeyType psrKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrainRate, nodeListKey);
  const KeyType meltKey = State<Dimension>::buildFieldKey(SolidFieldNames::meltSpecificEnergy, nodeListKey);
  CHECK(state.registered(PKey));
  CHECK(state.registered(epsKey));
  CHECK(state.registered(stressKey));
  CHECK(state.registered(psrKey));
  CHECK(state.registered(meltKey));

  const auto& P = state.field(PKey, 0.0);
  const auto& eps = state.field(epsKey, 0.0);
  const auto& S = state.field(stressKey, SymTensor::zero);
  const auto& psr = state.field(psrKey, 0.0);
  const auto& epsMelt = state.field(meltKey, 0.0);

  // Iterate over the internal nodes.
  const auto n = efail.numInternalElements();
#pragma omp parallel for
  for (auto i = 0; i < n; ++i) {
    if (-P(i) > -msigmamax) {
      efail(i) = mefailmin;
    } else {
      const auto sigmaVMinv = safeInv(sqrt(1.5*S(i).doubledot(S(i))));
      const auto chi = -P(i)*sigmaVMinv;
      if (chi > -mTcrit and chi < -msigmamax*sigmaVMinv) {
        const auto efailTcrit = (mD1(i) + mD2(i)*exp(-mD3*mTcrit))*
                                (1.0 + mD4*log(psr(i)*safeInv(mepsilondot0)))*
                                (1.0 + mD5*eps(i)*safeInv(epsMelt(i)));
        const auto psi = (chi + mTcrit)*safeInvVar(-msigmamax*sigmaVMinv + mTcrit);
        CHECK(psi >= 0.0 and psi <= 1.0);
        efail(i) = (1.0 - psi)*efailTcrit + psi*mefailmin;
      } else {
        efail(i) = (mD1(i) + mD2(i)*exp(-mD3*chi))*
                   (1.0 + mD4*log(psr(i)*safeInv(mepsilondot0)))*
                   (1.0 + mD5*eps(i)*safeInv(epsMelt(i)));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
JohnsonCookFailureStrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const JohnsonCookFailureStrainPolicy<Dimension>* rhsPtr = dynamic_cast<const JohnsonCookFailureStrainPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}
