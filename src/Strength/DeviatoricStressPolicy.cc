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
// Helper method to compute the J2 constant from the deviatoric stress.
//------------------------------------------------------------------------------
namespace {
template<typename Dimension>
inline
double
computeJ2(const typename Dimension::SymTensor& x) {
  return 0.5*x.doubledot(x);
}
}

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
  // FieldList<Dimension, Scalar> ps = state.fields(SolidFieldNames::plasticStrain, 0.0);
  // FieldList<Dimension, Scalar> psr = derivs.fields(SolidFieldNames::plasticStrainRate, 0.0);
  // FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  // const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  // const FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  // const FieldList<Dimension, Scalar> G = state.fields(SolidFieldNames::shearModulus, 0.0);
  // const FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  // const FieldList<Dimension, Scalar> ps0 = state.fields(SolidFieldNames::plasticStrain + "0", 0.0);
  // const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  const FieldList<Dimension, SymTensor> DSDt = derivs.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + 
                                                             SolidFieldNames::deviatoricStress, SymTensor::zero);

  // Iterate over the internal nodes.
  const unsigned numFields = S.numFields();
  for (unsigned k = 0; k != numFields; ++k) {
    const unsigned n = S[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {

      const SymTensor Sold = S(k,i);                        // Starting deviatoric stress.
      const SymTensor S0 = Sold + multiplier*(DSDt(k,i));    // Elastic prediction for the new deviatoric stress.

      // // Equivalent stress deviator.
      // const double J2 = computeJ2<Dimension>(S0);
      // CHECK(J2 >= 0.0);

      // // Radial return correction.
      // const double f = max(1.0, 3.0*J2*safeInvVar(FastMath::square(0.1*Y(k,i))));
      // CHECK(f >= 1.0);

      // // Check for yielding.
      // if (f > 1.0) {

      //   const SymTensor deformation = DvDx(k,i).Symmetric();
      //   const SymTensor linearDeformation = deformation - (deformation.Trace()/Dimension::nDim)*SymTensor::one;
      //   const Tensor spin = DvDx(k,i).SkewSymmetric();
      //   const SymTensor S1 = S0/sqrt(f);                                                       // New deviatoric stress with plastic yielding.
      //   const Tensor R = (spin*(Sold + S1)).SkewSymmetric();
      //   const Tensor deltaPStensor = linearDeformation*multiplier - 0.5*(S1 - Sold + R*multiplier)/G(k,i);
      //   const double deltaEPS = sqrt(2.0/3.0*deltaPStensor.doubledot(deltaPStensor));
      //   const double epsdot = sqrt(2.0/3.0*linearDeformation.doubledot(linearDeformation));
      //   const double plasticWork = Sold.doubledot(deltaPStensor);
      //   ps(k,i) += deltaEPS;
      //   psr(k,i) = (ps(k,i) - ps0(k,i))*safeInv(dt);
      //   // eps(k,i) += plasticWork/rho(k,i); // Should add plastic work to intenrnal energy?

      // } else {

        // Purely elastic flow.
        S(k,i) = S0;

      // }
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

