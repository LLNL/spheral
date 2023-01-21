#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// etamin_solid
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
etamin_solid() const {
  return mEtaMinSolid;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
etamin_solid(double x) {
  mEtaMinSolid = x;
}

//------------------------------------------------------------------------------
// etamax_solid
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
etamax_solid() const {
  return mEtaMaxSolid;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
etamax_solid(double x) {
  mEtaMaxSolid = x;
}

//------------------------------------------------------------------------------
// a
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
a() const {
  return ma;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
a(double x) {
  ma = x;
}

//------------------------------------------------------------------------------
// b
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
b() const {
  return mb;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
b(double x) {
  mb = x;
}

//------------------------------------------------------------------------------
// A
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
A(double x) {
  mA = x;
}

//------------------------------------------------------------------------------
// B
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
B(double x) {
  mB = x;
}

//------------------------------------------------------------------------------
// alpha
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
alpha() const {
  return malpha;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
alpha(double x) {
  malpha = x;
}

//------------------------------------------------------------------------------
// beta
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
beta() const {
  return mbeta;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
beta(double x) {
  mbeta = x;
}

//------------------------------------------------------------------------------
// eps0
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
eps0() const {
  return meps0;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
eps0(double x) {
  meps0 = x;
}

//------------------------------------------------------------------------------
// epsLiquid
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
epsLiquid() const {
  return mepsLiquid;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
epsLiquid(double x) {
  mepsLiquid = x;
}

//------------------------------------------------------------------------------
// epsVapor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
epsVapor() const {
  return mepsVapor;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
epsVapor(double x) {
  mepsVapor = x;
}

//------------------------------------------------------------------------------
// Access the atomic weight
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
TillotsonEquationOfState<Dimension>::
atomicWeight(double x) {
  mAtomicWeight = x;
}

//------------------------------------------------------------------------------
// The internal helper methods for computing components of the EOS.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
computePhi(const double& eta, const double& eps) const {
  return mb/(1.0 + eps/(meps0*eta*eta));
}

// P1
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
computeP1(const double& mu, const double& P2) const {
  return P2 + mB*mu*mu;
}

// P2
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
computeP2(const double& phi, const double& mu,
          const double& rho, const double& eps) const {
  return (ma + phi)*rho*eps + mA*mu;
}

// P4
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
computeP4(const double& phi, const double& mu, const double& eta,
          const double& rho, const double& eps) const {
  const double thpt = 1.0 - 1.0/eta;
  return ma*rho*eps + (phi*rho*eps + mA*mu*exp(mbeta*thpt))*exp(-malpha*thpt*thpt);
}

// dphi/drho|_eps
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dphidrho_eps(const double& rho0, const double& eta, const double& eps) const {
  return 2.0*mb*eps/(rho0*meps0*eta*eta*eta*FastMath::square(1.0 + eps/(meps0*eta*eta)));
}

// dP1/drho|_eps
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dP1drho_eps(const double& rho0, const double& mu, const double& dP2drho_eps) const {
  return dP2drho_eps + 2.0*mB*mu/rho0;
}

// dP2/drho|_eps
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dP2drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0,
                    const double& rho, const double& eps) const {
  return (ma + phi + rho*dphidrho_eps)*eps + mA/rho0;
}

// dP4/drho|_eps
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dP4drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0,
                    const double& eta, const double& mu,
                    const double& rho, const double& eps) const {
  const double thpt = 1.0 - 1.0/eta;
  return ma*eps +
    exp(-malpha*thpt*thpt)*((phi + rho*dphidrho_eps - 2.0*malpha*thpt*rho0/rho*phi)*eps +
                            mA*(mbeta*rho0/(rho*rho)*mu + 1.0/rho0 -
                                2.0*malpha*thpt*rho0/(rho*rho)*mu)*exp(mbeta*thpt));
//   return ma*eps +
//     exp(-malpha*thpt*thpt)*(2.0*malpha*thpt/rho0*(phi*rho*eps + mA*mu*exp(mbeta*thpt)) +
//                             phi*eps + rho*eps*dphidrho_eps +
//                             mA*(mbeta*rho0*mu*eta + 1.0/rho0)*exp(mbeta*thpt));
}

// dphi/deps|_rho
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dphideps_rho(const double& eta, const double& eps) const {
  return -mb/(meps0*eta*eta*FastMath::square(1.0 + eps/(meps0*eta*eta)));
}

// dP2/deps|_rho
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dP2deps_rho(const double& phi, const double& dphideps_rho,
                    const double& rho, const double& eps) const {
  return (ma + phi + eps*dphideps_rho)*rho;
}

// dP4/deps|_rho
template<typename Dimension>
inline
double
TillotsonEquationOfState<Dimension>::
compute_dP4deps_rho(const double& phi, const double& dphideps_rho, const double& eta,
                    const double& rho, const double& eps) const {
  const double thpt = 1.0 - 1.0/eta;
  return rho*(ma + (phi + eps*dphideps_rho)*exp(-malpha*thpt*thpt));
}

}
