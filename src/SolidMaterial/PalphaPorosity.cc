//---------------------------------Spheral++----------------------------------//
// PalphaPorosity -- Strain porosity EOS modifier.
// 
// See header for references and such.
//----------------------------------------------------------------------------//
#include "SolidMaterial/PalphaPorosity.hh"
#include "SolidMaterial/PalphaPressurePolicy.hh"
#include "FileIO/FileIO.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/IncrementFieldList.hh"

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
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
PalphaPorosity(PorousEquationOfState<Dimension>& porousEOS,
               PorousStrengthModel<Dimension>& porousStrength,
               const NodeList<Dimension>& nodeList,
               const double phi0,
               const double Pe,
               const double Pt,
               const double alphae,
               const double alphat,
               const double n1,
               const double n2,
               const double cS0,
               const double c0):
  Physics<Dimension>(),
  mPe(Pe),
  mPt(Pt),
  mAlphae(alphae),
  mAlphat(alphat),
  mn1(n1),
  mn2(n2),
  mcS0(cS0),
  mPorousEOS(porousEOS),
  mPorousStrength(porousStrength),
  mNodeList(nodeList),
  mc0(SolidFieldNames::porosityc0, nodeList, c0),
  mAlpha0(SolidFieldNames::porosityAlpha0, nodeList, 1.0/(1.0 - phi0)),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList, 1.0/(1.0 - phi0)),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mdPdU(HydroFieldNames::partialPpartialEps, nodeList),
  mdPdR(HydroFieldNames::partialPpartialRho, nodeList),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(mPe <= mPt,
          "PalphaPorosity input ERROR : require Pe <= Pt: Pe = " << mPe << ", Pt = " << mPt);
  const auto alpha0_min = mAlpha0.min();
  VERIFY2((1.0 <= mAlphae) and (mAlphat <= mAlphae) and (mAlphae <= alpha0_min),
          "PalphaPorosity input ERROR : require 1.0 <= alphat <= alphae <= alpha0, (alphat, alphae, alpha0) = " << mAlphat << ", " << mAlphae << ", " << alpha0_min);
  VERIFY2(phi0 >= 0.0 and phi0 < 1.0,
          "ERROR : Initial porosity required to be in the range phi0 = [0.0, 1.0) : phi0 = " << phi0);
  mPorousEOS.alpha(mAlpha);
  mPorousEOS.alpha0(mAlpha0);
  mPorousEOS.c0(mc0);
  mPorousStrength.alpha(mAlpha);
  ENSURE(mPorousEOS.valid());
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
PalphaPorosity(PorousEquationOfState<Dimension>& porousEOS,
               PorousStrengthModel<Dimension>& porousStrength,
               const NodeList<Dimension>& nodeList,
               const Field<Dimension, Scalar>& phi0,
               const double Pe,
               const double Pt,
               const double alphae,
               const double alphat,
               const double n1,
               const double n2,
               const double cS0,
               const Field<Dimension, Scalar>& c0):
  Physics<Dimension>(),
  mPe(Pe),
  mPt(Pt),
  mAlphae(alphae),
  mAlphat(alphat),
  mn1(n1),
  mn2(n2),
  mcS0(cS0),
  mPorousEOS(porousEOS),
  mPorousStrength(porousStrength),
  mNodeList(nodeList),
  mc0(SolidFieldNames::porosityc0, nodeList),
  mAlpha0(SolidFieldNames::porosityAlpha0, nodeList),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mdPdU(HydroFieldNames::partialPpartialEps),
  mdPdR(HydroFieldNames::partialPpartialRho),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(mPe <= mPt,
          "PalphaPorosity input ERROR : require Pe <= Pt: Pe = " << mPe << ", Pt = " << mPt);
  const auto phi0_min = phi0.min();
  const auto phi0_max = phi0.max();
  VERIFY2(phi0_min >= 0.0 and phi0_max < 1.0,
          "ERROR : Initial porosity required to be in the range phi0 = [0.0, 1.0): phi0 min/max = " << phi0_min << " " << phi0_max);
  const auto n = nodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    mc0[i] = c0[i];
    mAlpha0[i] = 1.0/(1.0 - phi0[i]);
    mAlpha[i] = 1.0/(1.0 - phi0[i]);
  }
  const auto alpha0_min = mAlpha0.min();
  VERIFY2((1.0 <= mAlphae) and (mAlphae <= mAlphat) and (mAlphat <= alpha0_min),
          "PalphaPorosity input ERROR : require 1.0 <= alphae <= alphat <= alpha0, (alphae, alphat, alpha0) = " << mAlphae << ", " << mAlphat << ", " << alpha0_min);
  mPorousEOS.alpha(mAlpha);
  mPorousEOS.alpha0(mAlpha0);
  mPorousEOS.c0(mc0);
  mPorousStrength.alpha(mAlpha);
  ENSURE(mPorousEOS.valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
~PalphaPorosity() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the alpha field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state fields.
  const auto  rhoKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, mNodeList.name());
  const auto  PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, mNodeList.name());
  const auto  dPduKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialEps, mNodeList.name());
  const auto  dPdrKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialRho, mNodeList.name());
  const auto  alphaKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, mNodeList.name());
  const auto  DalphaDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, mNodeList.name());
  const auto  DrhoDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + HydroFieldNames::massDensity, mNodeList.name());
  const auto  DuDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, mNodeList.name());
  const auto& rho = state.field(rhoKey, 0.0);
  const auto& P = state.field(PKey, 0.0);
  const auto& dPdu = state.field(dPduKey, 0.0);
  const auto& dPdr = state.field(dPdrKey, 0.0);
  const auto& alpha = state.field(alphaKey, 0.0);
  const auto& DrhoDt = derivs.field(DrhoDtKey, 0.0);
  const auto& DuDt = derivs.field(DuDtKey, 0.0);
  auto&       DalphaDt = derivs.field(DalphaDtKey, 0.0);

  const auto rho0 = mPorousEOS.referenceDensity();

  // Walk the nodes.
  const auto n = mNodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto rhoi = rho(i);
    const auto Pi = P(i);
    const auto dPsdui = dPdu(i);
    const auto dPsdri = dPdr(i);
    const auto alphai = alpha(i);
    const auto c0i = mc0(i);         // Initial sound speed with porosity
    const auto DrhoDti = DrhoDt(i);
    const auto DuDti = DuDt(i);
    const auto Psi = alphai*Pi;

    // Time evolution of the solid pressure
    const auto dPsdti = dPsdri*DrhoDti + dPsdui*DuDti;

    // First compute the derivative with respect to pressure
    Scalar DalphaDpi = 0.0;
    if (Pi < mPe or dPsdti < 0.0) {

      // Elastic
      if (c0i != mcS0) {  // If initial porous sound speed is the same as solid phase, no elastic evolution
        const auto halpha = 1.0 + (alphai - 1.0)*(c0i - mcS0)*safeInvVar(mcS0*(mAlphae - 1.0));
        DalphaDpi = alphai*alphai/(mcS0*mcS0*rho0)*(1.0 - safeInvVar(halpha*halpha));
      }

    } else {

      // Plastic
      DalphaDpi = (Pi < mPt ?
                   1.0 - mn1*(mAlphae - mAlphat)*pow((mPt - Pi)/(mPt - mPe), mn1)*safeInv(mPt - Pi) - mn2*(mAlphat - 1.0)*pow((Psi - Pi)/(Psi - mPe), mn2)*safeInv(Psi - Pi) :
                   1.0 - mn2*(mAlphat - 1.0)*pow((Psi - Pi)/(Psi - mPe), mn2)*safeInv(Psi - Pi));

    }

    // Now we can compute the final time derivative
    const auto dPdti = (alphai*dPsdri*DrhoDti + dPsdui*DuDti)*safeInvVar(alphai + DalphaDpi*(Pi - rhoi*dPsdri));
    DalphaDt(i) = DalphaDpi*dPdti;
  }
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename PalphaPorosity<Dimension>::TimeStepType
PalphaPorosity<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return TimeStepType(1.0e100, "Rate of porosity change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
registerState(DataBase<Dimension>& /*dataBase*/,
              State<Dimension>& state) {

  // Override the pressure policy to compute the partial derivatives of the 
  // pressure as well
  auto pressurePolicy = std::make_shared<PalphaPressurePolicy<Dimension>>();
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  state.enroll(pressure, pressurePolicy);
  state.enroll(mdPdU);
  state.enroll(mdPdR);

  // Register the P-alpha state
  auto alphaPolicy = std::make_shared<IncrementBoundedState<Dimension, Scalar, Scalar>>(1.0);
  state.enroll(mAlpha, alphaPolicy);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDalphaDt);
}

//------------------------------------------------------------------------------
// One time initializations at problem set up.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& /*dataBase*/) {
  // Initialize the distention field.
  mAlpha = mAlpha0;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mc0, pathName + "/c0");
  file.write(mAlpha0, pathName + "/alpha0");
  file.write(mAlpha, pathName + "/alpha");
  file.write(mDalphaDt, pathName + "/DalphaDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mc0, pathName + "/c0");
  file.read(mAlpha0, pathName + "/alpha0");
  file.read(mAlpha, pathName + "/alpha");
  file.read(mDalphaDt, pathName + "/DalphaDt");
}

//------------------------------------------------------------------------------
// Compute the current porosity
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
PalphaPorosity<Dimension>::
phi() const {
  Field<Dimension, Scalar> phi("porosity", mNodeList, 0.0);
  const auto n = mNodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    CHECK(mAlpha(i) > 0.0);
    phi(i) = 1.0 - 1.0*safeInvVar(mAlpha(i));
  }
  return phi;
}

}