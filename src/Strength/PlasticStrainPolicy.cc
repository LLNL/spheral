//---------------------------------Spheral++----------------------------------//
// PlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
//
// Created by JMO, Wed Sep 15 23:07:21 PDT 2004
//----------------------------------------------------------------------------//
#include "PlasticStrainPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"

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
namespace {

//------------------------------------------------------------------------------
// Helper method to compute the J2 constant from the deviatoric stress.
//------------------------------------------------------------------------------
inline
double
computeJ2(const Dim<3>::SymTensor& S) {
  return 0.5*(S.doubledot(S));
}

inline
double
computeJ2(const Dim<2>::SymTensor& S) {
  // the third diagonal component of S (acutally negative to make Tr(S)=0)
  const auto S33 = S.Trace();
  return 0.5*(S.doubledot(S) + S33*S33);
}

inline
double
computeJ2(const Dim<1>::SymTensor& S) {
  // S_22 == S_33 = -S_11/2
  return 0.75*S.xx()*S.xx();
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PlasticStrainPolicy<Dimension>::
PlasticStrainPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(SolidFieldNames::deviatoricStress,
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
       const double /*multiplier*/,
       const double /*t*/,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::plasticStrain and
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto ps = state.fields(fieldKey, 0.0);
  const auto numFields = ps.numFields();

  // Get the state we depend on.
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto P = state.fields(HydroFieldNames::pressure, 0.0);
  const auto G = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  const auto ps0 = state.fields(SolidFieldNames::plasticStrain + "0", 0.0);
  auto       deviatoricStress = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto       psr = derivs.fields(SolidFieldNames::plasticStrainRate, 0.0);

  // Iterate over the internal nodes.
  for (auto k = 0u; k < numFields; ++k) {
    const auto n = ps[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {

      // Equivalent stress deviator.
      const auto J2 = computeJ2(deviatoricStress(k,i));
      CHECK(J2 >= 0.0);
      const auto equivalentStressDeviator = sqrt(3.0*J2);

      // Plastic yield limit.
      const auto yieldLimit = max(0.0, Y(k,i));

      // von Mises yield scaling constant.
      const auto f = min(1.0, yieldLimit*safeInvVar(equivalentStressDeviator));

      // Scale the stress deviator.
      deviatoricStress(k,i) *= f;

      // Update the plastic strain and strain rate.
      if (distinctlyGreaterThan(G(k,i), 0.0)) {
        ps(k,i) += (1.0 - f)*equivalentStressDeviator / (3.0*G(k,i));
        psr(k,i) = (ps(k,i) - ps0(k,i))*safeInv(dt);
      }
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

