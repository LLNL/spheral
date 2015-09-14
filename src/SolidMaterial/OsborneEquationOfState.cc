//---------------------------------Spheral++----------------------------------//
// OsborneEquationOfState -- Osborne  equation of state.
//
// Reference: PAGOSA Physics manual, LA-14425-M
//----------------------------------------------------------------------------//

#include "OsborneEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with the given Osborne constants.
//------------------------------------------------------------------------------
template<typename Dimension>
OsborneEquationOfState<Dimension>::
OsborneEquationOfState(const double referenceDensity,
                       const double etamin,
                       const double etamax,
                       const double a1,
                       const double a2pos,
                       const double a2neg,
                       const double b0,
                       const double b1,
                       const double b2pos,
                       const double b2neg,
                       const double c0,
                       const double c1,
                       const double c2pos,
                       const double c2neg,
                       const double E0,
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
  mA1(a1),
  mA2pos(a2pos),
  mA2neg(a2neg),
  mB0(b0),
  mB1(b1),
  mB2pos(b2pos),
  mB2neg(b2neg),
  mC0(c0),
  mC1(c1),
  mC2pos(c2pos),
  mC2neg(c2neg),
  mE0(E0),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * constants.molarGasConstant() / mAtomicWeight),
  mExternalPressure(externalPressure) {
  VERIFY(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
OsborneEquationOfState<Dimension>::~OsborneEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  const double rho0 = this->referenceDensity();
  for (int i = 0; i != pressure.size(); ++i) {
    const double eta = this->boundedEta(massDensity(i));
    const double mu = eta - 1.0;
    const double E = rho0*specificThermalEnergy(i);
    const double a2 = mu > 0.0 ? mA2pos : mA2neg;
    const double b2 = mu > 0.0 ? mB2pos : mB2neg;
    const double c2 = mu > 0.0 ? mC2pos : mC2neg;
    pressure(i) = this->applyPressureLimits((mA1*mu + a2*mu*mu +
                                             (mB0 + mB1*mu + b2*mu*mu)*E +
                                             (mC0 + mC1*mu + c2*mu*mu)*E*E)*safeInvVar(E + mE0) -
                                            mExternalPressure);
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  VERIFY2(false, "temperature unimplemented for the Osborne equation of state.");
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  VERIFY2(false, "specific thermal energy unimplemented for the Osborne equation of state.");
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (int i = 0; i != soundSpeed.size(); ++i) {
    const double dPdrhoi = DPDrho(massDensity(i), specificThermalEnergy(i));
    CHECK(dPdrhoi >= 0.0);
    soundSpeed(i) = sqrt(dPdrhoi);
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  VERIFY2(false, "gamma unimplemented for the Osborne equation of state.");
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
OsborneEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (int i = 0; i != bulkModulus.size(); ++i) {
    bulkModulus(i) = massDensity(i) * DPDrho(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho).
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension>
double
OsborneEquationOfState<Dimension>::
DPDrho(const double massDensity,
       const double specificThermalEnergy) const {
  const double tiny = 1.0e-20;
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double E = rho0*specificThermalEnergy;
  const double a2 = mu > 0.0 ? mA2pos : mA2neg;
  const double b2 = mu > 0.0 ? mB2pos : mB2neg;
  const double c2 = mu > 0.0 ? mC2pos : mC2neg;
  const double P = (mA1*mu + a2*mu*mu +
                    (mB0 + mB1*mu + b2*mu*mu)*E +
                    (mC0 + mC1*mu + c2*mu*mu)*E*E)*safeInvVar(E + mE0);
  const double DPDrho_e = (mA1 + 2.0*a2*mu + E*(mB1 + 2.0*b2*mu + mC1*E))*safeInvVar(rho0*(E + mE0));
  const double DPDeps_r = rho0*(mB0 + mB1*mu + b2*mu*mu + 2.0*E*(mC0 + mC1*mu) - P)*safeInv(E + mE0);
  return max(0.0, DPDrho_e + P/(rho*rho)*DPDeps_r);
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
OsborneEquationOfState<Dimension>::valid() const {
  return true;
}

}
}

