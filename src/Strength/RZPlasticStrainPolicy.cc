//---------------------------------Spheral++----------------------------------//
// RZPlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
// This one is also specialized for the RZ strength case.
//
// Created by JMO, Mon May  9 14:09:12 PDT 2016
//----------------------------------------------------------------------------//
#include "RZPlasticStrainPolicy.hh"
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

//------------------------------------------------------------------------------
// Helper method to compute the J2 constant from the deviatoric stress.
//------------------------------------------------------------------------------
inline
double
computeJ2(const Dim<2>::SymTensor& x,
          const Dim<2>::Scalar& xTT) {
  return 0.5*(x.doubledot(x) + xTT*xTT);
}


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RZPlasticStrainPolicy::
RZPlasticStrainPolicy():
  FieldListUpdatePolicyBase<Dim<2>, Dim<2>::Scalar>(SolidFieldNames::deviatoricStress,
                                                    HydroFieldNames::massDensity,
                                                    HydroFieldNames::specificThermalEnergy,
                                                    HydroFieldNames::pressure) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
RZPlasticStrainPolicy::
~RZPlasticStrainPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
// Note this updates both the deviatoric stress and the plastic strain
// via the von Mises yielding correction.
//------------------------------------------------------------------------------
void
RZPlasticStrainPolicy::
update(const KeyType& key,
       State<Dim<2> >& state,
       StateDerivatives<Dim<2> >& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::plasticStrain and
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> ps = state.fields(fieldKey, 0.0);
  const unsigned numFields = ps.numFields();

  // Get the state we depend on.
  auto       deviatoricStress = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto P = state.fields(HydroFieldNames::pressure, 0.0);
  const auto G = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  const auto ps0 = state.fields(SolidFieldNames::plasticStrain + "0", 0.0);
  auto       psr = derivs.fields(SolidFieldNames::plasticStrainRate, 0.0);

  // Iterate over the internal nodes.
  for (unsigned k = 0; k != numFields; ++k) {
    const unsigned n = ps[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {

      // Equivalent stress deviator.
      const double J2 = computeJ2(deviatoricStress(k,i), -deviatoricStress(k,i).Trace());
      CHECK(J2 >= 0.0);
      const double equivalentStressDeviator = sqrt(3.0*J2);

      // Plastic yield limit.
      const double yieldLimit = max(0.0, Y(k,i));

      // von Mises yield scaling constant.
      double f;
      if (distinctlyGreaterThan(equivalentStressDeviator, 0.0)) {
        f = min(1.0, yieldLimit/equivalentStressDeviator);
      } else {
        f = 1.0;
      }

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
bool
RZPlasticStrainPolicy::
operator==(const UpdatePolicyBase<Dim<2> >& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const RZPlasticStrainPolicy* rhsPtr = dynamic_cast<const RZPlasticStrainPolicy*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

