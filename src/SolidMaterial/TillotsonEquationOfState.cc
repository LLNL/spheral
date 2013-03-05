//---------------------------------Spheral++----------------------------------//
// TillotsonEquationOfState -- This equation of state is designed to represent
// metallic materials over a range of pressure and density -- spanning solid, 
// liquid, and vapor states.
//
// Tillotson 1962
//
// Created by JMO, Wed Mar 16 23:31:17 PDT 2011
//----------------------------------------------------------------------------//
#include <iostream>
using namespace std;

#include "TillotsonEquationOfState.hh"
#include "Field/Field.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with the given coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TillotsonEquationOfState<Dimension>::
TillotsonEquationOfState(const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double a,
                         const double b,
                         const double A,
                         const double B,
                         const double alpha,
                         const double beta,
                         const double eps0,
                         const double epsLiquid,
                         const double epsVapor,
                         const double atomicWeight,
                         const Material::PhysicalConstants& constants,
                         const double externalPressure,
                         const double minimumPressure,
                         const double maximumPressure):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure),
  ma(a),
  mb(b),
  mA(A),
  mB(B),
  malpha(alpha),
  mbeta(beta),
  meps0(eps0),
  mepsLiquid(epsLiquid),
  mepsVapor(epsVapor),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * referenceDensity * constants.molarGasConstant() / atomicWeight),
  mExternalPressure(externalPressure) {
  VERIFY(distinctlyGreaterThan(mAtomicWeight, 0.0));
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TillotsonEquationOfState<Dimension>::
~TillotsonEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != Pressure.size(); ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != temperature.size(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  for (int i = 0; i != specificThermalEnergy.size(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != soundSpeed.size(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  VERIFY2(false, "gamma not defined for Tillotson EOS!");
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TillotsonEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != bulkModulus.size(); ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.

  // Define three fundamental pressures:
  //   P1 - solid, compression.
  //   P2 - solid, expansion.
  //   P4 - gaseous, expansion.
  // A third category (P3) is interpreted as a mixture of gaseous and solid
  // phases, and is interpolated between P2 and P3.

  // There are four regimes:
  double P;
  if (mu >= 0.0) {

    // Option 1: compression, solid.
    const double phi = computePhi(eta, eps);
    P = computeP1(mu, computeP2(phi, mu, rho, eps));

  } else if (eps <= mepsLiquid) {

    // Option 2: expansion, solid : same as 1, but setting B=0.
    const double phi = computePhi(eta, eps);
    P = computeP2(phi, mu, rho, eps);

  } else if (eps <= mepsVapor) {

    // Option 3: expansion, liquid.
    // Treated here as a linear combination of the solid and gaseous phases.
    // Following Saito et al. we compute P2 and P4 at the epsLiquid and epsVapor 
    // specific energies.
    const double phi2 = computePhi(eta, mepsLiquid);
    const double phi4 = computePhi(eta, mepsVapor);
    const double P2 = computeP2(phi2, mu, rho, mepsLiquid);
    const double P4 = computeP4(phi4, mu, eta, rho, mepsVapor);
    P = P2 + (P4 - P2)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);

  } else {

    // Option 4: expansion, gaseous.
    const double phi = computePhi(eta, eps);
    P = computeP4(phi, mu, eta, rho, eps);

  }

  // That's it.
  return max(this->minimumPressure(), min(this->maximumPressure(), P - mExternalPressure));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  const double eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.
  return mCv*eps + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  return (temperature - 300.0)/mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  const double c2 = computeDPDrho(massDensity, specificThermalEnergy);
  //ENSURE(c2 >= 0.0);
  return sqrt(max(0.0, c2));
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  VERIFY2(false, "gamma not defined for Tillotson EOS!");
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho) (adiabatic)
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension>
double
TillotsonEquationOfState<Dimension>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.
  const double Prho2 = this->pressure(massDensity, specificThermalEnergy)/(rho*rho);

  // There are four regimes:
  double dPdrho;
  if (mu >= 0.0) {

    // Option 1: compression, solid.
    const double phi = computePhi(eta, eps);
    const double dphidrho_eps = compute_dphidrho_eps(rho0, eta, eps);
    const double dphideps_rho = compute_dphideps_rho(eta, eps);
    const double dP2drho_eps = compute_dP2drho_eps(phi, dphidrho_eps, rho0, rho, eps);
    const double dP1drho_eps = compute_dP1drho_eps(rho0, mu, dP2drho_eps);
    const double dP2deps_rho = compute_dP2deps_rho(phi, dphideps_rho, rho, eps);
    dPdrho = dP1drho_eps + Prho2*dP2deps_rho; // Note dP1deps_rho == dP2deps_rho

  } else if (eps <= mepsLiquid) {

    // Option 2: expansion, solid : same as 1, but setting B=0.
    const double phi = computePhi(eta, eps);
    const double dphidrho_eps = compute_dphidrho_eps(rho0, eta, eps);
    const double dphideps_rho = compute_dphideps_rho(eta, eps);
    const double dP2drho_eps = compute_dP2drho_eps(phi, dphidrho_eps, rho0, rho, eps);
    const double dP2deps_rho = compute_dP2deps_rho(phi, dphideps_rho, rho, eps);
    dPdrho = dP2drho_eps + Prho2*dP2deps_rho;

  } else if (eps <= mepsVapor) {

    // Option 3: expansion, liquid.
    // Treated here as a linear combination of the solid and gaseous phases.
    // Following Saito et al. we compute P2 and P4 at the epsLiquid and epsVapor 
    // specific energies.
    double dP2drho, dP4drho;
    {
      const double phi = computePhi(eta, mepsLiquid);
      const double dphidrho_eps = compute_dphidrho_eps(rho0, eta, mepsLiquid);
      const double dphideps_rho = compute_dphideps_rho(eta, mepsLiquid);
      const double dP2drho_eps = compute_dP2drho_eps(phi, dphidrho_eps, rho0, rho, mepsLiquid);
      const double dP2deps_rho = compute_dP2deps_rho(phi, dphideps_rho, rho, mepsLiquid);
      dP2drho = dP2drho_eps + Prho2*dP2deps_rho;
    }
    {
      const double phi = computePhi(eta, mepsVapor);
      const double dphidrho_eps = compute_dphidrho_eps(rho0, eta, mepsVapor);
      const double dphideps_rho = compute_dphideps_rho(eta, mepsVapor);
      const double dP4drho_eps = compute_dP4drho_eps(phi, dphidrho_eps, rho0, eta, mu, rho, mepsVapor);
      const double dP4deps_rho = compute_dP4deps_rho(phi, dphideps_rho, eta, rho, mepsVapor);
      dP4drho = dP4drho_eps + Prho2*dP4deps_rho;
    }
    dPdrho = dP2drho + (dP4drho - dP2drho)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);

  } else {

    // Option 4: expansion, gaseous.
    const double phi = computePhi(eta, eps);
    const double dphidrho_eps = compute_dphidrho_eps(rho0, eta, eps);
    const double dphideps_rho = compute_dphideps_rho(eta, eps);
    const double dP4drho_eps = compute_dP4drho_eps(phi, dphidrho_eps, rho0, eta, mu, rho, eps);
    const double dP4deps_rho = compute_dP4deps_rho(phi, dphideps_rho, eta, rho, eps);
    dPdrho = dP4drho_eps + Prho2*dP4deps_rho;

  }

  // That's it.
  return std::abs(dPdrho);
}

}
}

