//---------------------------------Spheral++----------------------------------//
// PlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
//
// Created by JMO, Wed Sep 15 23:07:21 PDT 2004
//----------------------------------------------------------------------------//
#include "PlasticStrainPolicy.hh"
#include "SolidNodeList.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

using namespace std;

using NodeSpace::NodeList;
using FieldSpace::Field;
using SolidMaterial::SolidNodeList;

//------------------------------------------------------------------------------
// Helper method to compute the J2 constant from the deviatoric stress.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
computeJ2(const typename Dimension::SymTensor& x) {
  return 0.5*x.doubledot(x);
}

// template<>
// inline
// double
// computeJ2< Dim<1> >(const Dim<1>::SymTensor& x) {
//   const double x22_33 = 0.5*(x.xx());
//   return 0.5*(x.doubledot(x) + 2.0*x22_33*x22_33);
// }

// template<>
// inline
// double
// computeJ2< Dim<2> >(const Dim<2>::SymTensor& x) {
//   const double x33 = x.Trace();
//   return 0.5*(x.doubledot(x) + x33*x33);
// }

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PlasticStrainPolicy<Dimension>::
PlasticStrainPolicy():
  UpdatePolicyBase<Dimension>(SolidFieldNames::deviatoricStress,
                              HydroFieldNames::massDensity,
                              HydroFieldNames::specificThermalEnergy,
                              HydroFieldNames::pressure) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PlasticStrainPolicy<Dimension>::
~PlasticStrainPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
// Note this updates both the deviatoric stress and the plastic strain
// via the von Mises yielding correction.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlasticStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::plasticStrain);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // Get the deviatoric stress.
  const KeyType deviatoricStressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  CHECK(state.registered(deviatoricStressKey));
  Field<Dimension, SymTensor>& deviatoricStress = state.field(deviatoricStressKey, SymTensor::zero);

  // We also need the mass density, thermal energy, pressure, and plastic strain rate.
  const KeyType rhoKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType epsKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType GKey = State<Dimension>::buildFieldKey(SolidFieldNames::shearModulus, nodeListKey);
  const KeyType YKey = State<Dimension>::buildFieldKey(SolidFieldNames::yieldStrength, nodeListKey);
  const KeyType ps0Key = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrain + "0", nodeListKey);
  const KeyType psrKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrainRate, nodeListKey);
  CHECK(state.registered(rhoKey));
  CHECK(state.registered(epsKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(GKey));
  CHECK(state.registered(YKey));
  CHECK(state.registered(ps0Key));
  CHECK(derivs.registered(psrKey));
  const Field<Dimension, Scalar>& rho = state.field(rhoKey, 0.0);
  const Field<Dimension, Scalar>& eps = state.field(epsKey, 0.0);
  const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);
  const Field<Dimension, Scalar>& G = state.field(GKey, 0.0);
  const Field<Dimension, Scalar>& Y = state.field(YKey, 0.0);
  const Field<Dimension, Scalar>& ps0 = state.field(ps0Key, 0.0);
  Field<Dimension, Scalar>& psr = derivs.field(psrKey, 0.0);

  // Iterate over the internal nodes.
  const SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateField.nodeListPtr());
  CHECK(solidNodeListPtr != 0);
  for (int i = 0; i != solidNodeListPtr->numInternalNodes(); ++i) {

    // Equivalent stress deviator.
    const double J2 = computeJ2<Dimension>(deviatoricStress(i));
    CHECK(J2 >= 0.0);
    const double equivalentStressDeviator = sqrt(3.0*J2);

    // Plastic yield limit.
    const double yieldLimit = max(0.0, Y(i));

    // von Mises yield scaling constant.
    double f;
    if (distinctlyGreaterThan(equivalentStressDeviator, 0.0)) {
      f = min(1.0, yieldLimit/equivalentStressDeviator);
    } else {
      f = 1.0;
    }

    // Scale the stress deviator.
    deviatoricStress(i) *= f;

    // Update the plastic strain and strain rate.
    if (distinctlyGreaterThan(G(i), 0.0)) {
      stateField(i) += (1.0 - f)*equivalentStressDeviator / (3.0*G(i));
      psr(i) = (stateField(i) - ps0(i))*safeInv(dt);
    }
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PlasticStrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const PlasticStrainPolicy<Dimension>* rhsPtr = dynamic_cast<const PlasticStrainPolicy<Dimension>*>(&rhs);
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
  template class PlasticStrainPolicy<Dim<1> >;
  template class PlasticStrainPolicy<Dim<2> >;
  template class PlasticStrainPolicy<Dim<3> >;
}

