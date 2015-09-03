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
                         const double etamin_solid,
                         const double etamax_solid,
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
                         const double maximumPressure,
                         const Material::MaterialPressureMinType minPressureType):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  mEtaMinSolid(etamin_solid),
  mEtaMaxSolid(etamax_solid),
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
  mCv(3.0 * constants.molarGasConstant() / atomicWeight),
  mExternalPressure(externalPressure) {
  VERIFY(mAtomicWeight>0.0);
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
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
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
  const double eta = this->boundedEta(massDensity),
               mu = eta - 1.0,
               rho0 = this->referenceDensity(),
               rho = rho0*eta,
               eps = std::max(0.0, specificThermalEnergy);   // I'm not sure if this EOS admits negative energies.

  // Tillotson defines three fundamental pressure regimes:
  //   P1 - solid, compression.
  //   P2 - solid, expansion.
  //   P4 - gaseous, expansion.
  // A forth category (P3) is interpreted as a mixture of gaseous and solid
  // phases, and is interpolated between P2 and P4.

  double P;
  const double phi = mb/(1.0 + eps/(meps0*eta*eta));
  const double chi = 1.0/eta - 1.0;

  if (mu >= 0.0) {

    // Regime 1: compression, solid.
    P = (ma + phi)*rho*eps + mA*mu + mB*mu*mu;

  } else if (eps <= mepsLiquid) {

    // Regime 2: expansion, solid : same as 1, only if rho>cutoff density.
    P = (eta > mEtaMinSolid) ? (ma + phi)*rho*eps + mA*mu + mB*mu*mu : 0.0;

  } else if (eps >= mepsVapor) {

    // Regime 4: expansion, vapor.
    P = ma*rho*eps + (phi*rho*eps + mA*mu*exp(-mbeta*chi))*exp(-malpha*chi*chi);

  } else {

    // Regime 3: expansion, partial vapor.
    // Treated here as a linear combination of the solid and gaseous phases.
    // Following <strike>Saito et al. we compute P2 and P4 at the epsLiquid and
    // epsVapor specific energies</strike> Melosh (personal communication) we
    // compute P2 and P4 at the given energy.
    double P2 = (eta > mEtaMinSolid) ? (ma + phi)*rho*eps + mA*mu + mB*mu*mu : 0.;
    double P4 = ma*rho*eps + 
               (phi*rho*eps + mA*mu*exp(-mbeta*chi))*exp(-malpha*chi*chi);
    P = P2 + (P4 - P2)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);

  }

  // That's it.
  return this->applyPressureLimits(P - mExternalPressure);
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
  return eps/mCv + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TillotsonEquationOfState<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  return (temperature - 300.0)*mCv;
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
  const double eta = this->boundedEta(massDensity),
               mu = eta - 1.0,
               rho0 = this->referenceDensity(),
               rho = rho0*eta,
               eps = std::max(0.0, specificThermalEnergy);

  double dPdrho_ad, dPdrho_eps, dPdeps_rho;
  const double phi = mb/(1.0 + eps/(meps0*eta*eta));
  const double chi = 1./eta - 1.0;
  const double dphidrho = (1./rho0)*((2.0*mb*meps0*eps*eta)/
                                     ((eps+meps0*eta*eta)*(eps+meps0*eta*eta)));
  const double dphideps = -mb*meps0*eta*eta/((eps+meps0*eta*eta)*(eps+meps0*eta*eta));

  if (mu >= 0.0) {

    // Regime 1: compression, solid.
    dPdrho_eps = eps*(ma + phi + rho*dphidrho) + (1./rho0)*(mA + 2.0*mB*mu);
    dPdeps_rho = rho*(ma + phi + eps*dphideps);

  } else if (eps <= mepsLiquid) {

    // Regime 2: expansion, solid : same as 1, but only if rho>cutoff density
    dPdrho_eps = (eta > mEtaMinSolid) ? 
              eps*(ma + phi + rho*dphidrho) + (1./rho0)*(mA + 2.0*mB*mu) : 0.0;
    dPdeps_rho = (eta > mEtaMinSolid) ?
              rho*(ma + phi + eps*dphideps) : 0.0;

  } else if (eps >= mepsVapor) {
     
    // Regime 4: expansion, gaseous.
    dPdrho_eps = ma*eps + 
             eps*exp(-malpha*chi*chi)*(phi + rho*dphidrho + 2.0*malpha*chi*phi*rho0/rho) + 
             mA*exp(-(malpha*chi*chi+mbeta*chi))*(1./rho0 + rho0*mu*mbeta/(rho*rho) + 
                                                        2.0*rho0*mu*malpha*chi/(rho*rho));
    dPdeps_rho = rho*(ma + exp(-malpha*chi*chi)*(phi + eps*dphideps));

  } else {

    // Regime 3: expansion, liquid.
    // Treated here as a linear combination of the solid and gaseous phases.
    double dP2drho_eps, dP4drho_eps, dP2deps_rho, dP4deps_rho;

    dP2drho_eps = (eta > mEtaMinSolid) ? 
               eps*(ma + phi + rho*dphidrho) + (1./rho0)*(mA + 2.0*mB*mu) : 0.0;
    dP4drho_eps = ma*eps + 
              eps*exp(-malpha*chi*chi)*(phi + rho*dphidrho + 2.0*malpha*chi*phi*rho0/rho) + 
              mA*exp(-(malpha*chi*chi+mbeta*chi))*(1./rho0 + rho0*mu*mbeta/(rho*rho) + 
                                                         2.0*rho0*mu*malpha*chi/(rho*rho));
    dP2deps_rho = (eta > mEtaMinSolid) ?
               rho*(ma + phi + eps*dphideps) : 0.0;
    dP4deps_rho = rho*(ma + exp(-malpha*chi*chi)*(phi + eps*dphideps));

    dPdrho_eps = dP2drho_eps + (dP4drho_eps - dP2drho_eps)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);
    dPdeps_rho = dP2deps_rho + (dP4deps_rho - dP2deps_rho)*(eps - mepsLiquid)/(mepsVapor - mepsLiquid);

  }

  // That's it.
  dPdrho_ad = dPdrho_eps + dPdeps_rho*pressure(rho, eps)/(rho*rho);
  return std::abs(dPdrho_ad);
}

}
}

