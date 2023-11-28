//---------------------------------Spheral++----------------------------------//
// PorosityModel
// Base class for PorosityModels for common functionality.
//
// Created by JMO, Mon Nov 13 09:28:26 PST 2023
//----------------------------------------------------------------------------//
#include "Porosity/PorosityModel.hh"
#include "Porosity/PorositySolidMassDensityPolicy.hh"
#include "Porosity/PorousBulkModulusPolicy.hh"
#include "Porosity/PorousSoundSpeedPolicy.hh"
#include "Porosity/PorousGammaPolicy.hh"
#include "Porosity/PorousEntropyPolicy.hh"
#include "Porosity/PorousShearModulusPolicy.hh"
#include "Porosity/PorousYieldStrengthPolicy.hh"
#include "Material/EquationOfState.hh"
#include "FileIO/FileIO.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"

#include <string>

using std::vector;
using std::string;
using std::to_string;
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
PorosityModel<Dimension>::
PorosityModel(const SolidNodeList<Dimension>& nodeList,
              const double phi0,
              const double cS0,
              const double c0,
              const double rhoS0,
              const bool jutziStateUpdate):
  Physics<Dimension>(),
  mJutziStateUpdate(jutziStateUpdate),
  mRhoS0(rhoS0),
  mcS0(cS0),
  mKS0(rhoS0*cS0*cS0),
  mfdt(0.5),
  mMaxAbsDalphaDt(0.0),
  mNodeList(nodeList),
  mAlpha0(SolidFieldNames::porosityAlpha0, nodeList, 1.0/(1.0 - phi0)),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList, 1.0/(1.0 - phi0)),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mSolidMassDensity(SolidFieldNames::porositySolidDensity, nodeList),
  mc0(SolidFieldNames::porosityc0, nodeList, c0),
  mfDSptr(),
  mfDSnewPtr(),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(phi0 >= 0.0 and phi0 < 1.0,
          "ERROR : Initial porosity required to be in the range phi0 = [0.0, 1.0) : phi0 = " << phi0);
  if (mJutziStateUpdate) {
    mfDSptr = std::make_shared<Field<Dimension, Scalar>>(SolidFieldNames::fDSjutzi, nodeList, 1.0);
    mfDSnewPtr = std::make_shared<Field<Dimension, Scalar>>(ReplaceBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::fDSjutzi, nodeList, 1.0);
  }
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorosityModel<Dimension>::
PorosityModel(const SolidNodeList<Dimension>& nodeList,
              const Field<Dimension, Scalar>& phi0,
              const double cS0,
              const Field<Dimension, Scalar>& c0,
              const double rhoS0,
              const bool jutziStateUpdate):
  Physics<Dimension>(),
  mJutziStateUpdate(jutziStateUpdate),
  mRhoS0(rhoS0),
  mcS0(cS0),
  mKS0(rhoS0*cS0*cS0),
  mfdt(0.5),
  mMaxAbsDalphaDt(0.0),
  mNodeList(nodeList),
  mAlpha0(SolidFieldNames::porosityAlpha0, nodeList),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mSolidMassDensity(SolidFieldNames::porositySolidDensity, nodeList),
  mc0(c0),
  mfDSptr(),
  mfDSnewPtr(),
  mRestart(registerWithRestart(*this)) {
  const auto phi0_min = phi0.min();
  const auto phi0_max = phi0.max();
  VERIFY2(phi0_min >= 0.0 and phi0_max < 1.0,
          "ERROR : Initial porosity required to be in the range phi0 = [0.0, 1.0): phi0 min/max = " << phi0_min << " " << phi0_max);
  const auto n = nodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    mAlpha0[i] = 1.0/(1.0 - phi0[i]);
    mAlpha[i] = 1.0/(1.0 - phi0[i]);
  }
  if (mJutziStateUpdate) {
    mfDSptr = std::make_shared<Field<Dimension, Scalar>>(SolidFieldNames::fDSjutzi, nodeList, 1.0);
    mfDSnewPtr = std::make_shared<Field<Dimension, Scalar>>(ReplaceBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::fDSjutzi, nodeList, 1.0);
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorosityModel<Dimension>::
~PorosityModel() {
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename PorosityModel<Dimension>::TimeStepType
PorosityModel<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

  // Just limit by fractional change (assume we're comparing with min(alpha) = 1.0)
  const auto dt = (mfdt > 0.0 ?
                   mfdt * safeInvVar(mMaxAbsDalphaDt) :
                   std::numeric_limits<double>::max());
  return TimeStepType(dt,
                      "Rate of porosity change: max(DalphaDt) = " + to_string(mMaxAbsDalphaDt) + "\n" +
                      "                              material = " + mNodeList.name() + "\n" +
                      "                               on rank = " + to_string(Process::getRank()));
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorosityModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Alias for shorter call building State Field keys
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, mNodeList.name()); };

  // Register the solid mass density
  state.enroll(mSolidMassDensity, make_policy<PorositySolidMassDensityPolicy<Dimension>>());

  // Register the distension
  state.enroll(mAlpha, (mJutziStateUpdate ?
                        make_policy<IncrementBoundedState<Dimension, Scalar, Scalar>>({SolidFieldNames::deviatoricStress, SolidFieldNames::fDSjutzi}, 1.0) :
                        make_policy<IncrementBoundedState<Dimension, Scalar, Scalar>>(1.0)));
  state.enroll(mAlpha0);

  // Initial sound speed
  state.enroll(mc0);

  // Check what other state is registered which needs to be overridden for porosity
  auto optionalOverridePolicy = [&](const std::string& fkey, std::shared_ptr<UpdatePolicyBase<Dimension>> policy) -> void {
                                  const auto fullkey = buildKey(fkey);
                                  if (state.registered(fullkey)) {
                                    auto& f = state.field(fullkey, 0.0);
                                    state.enroll(f, policy);
                                  }
                                };
  // optionalOverridePolicy(SolidFieldNames::bulkModulus, std::make_shared<PorousBulkModulusPolicy<Dimension>>());
  // optionalOverridePolicy(HydroFieldNames::soundSpeed, std::make_shared<PorousSoundSpeedPolicy<Dimension>>());
  optionalOverridePolicy(HydroFieldNames::gamma, make_policy<PorousGammaPolicy<Dimension>>());
  optionalOverridePolicy(HydroFieldNames::entropy, make_policy<PorousEntropyPolicy<Dimension>>());
  optionalOverridePolicy(SolidFieldNames::shearModulus, make_policy<PorousShearModulusPolicy<Dimension>>());

  // The state update descibed in Jutzi et al. 2008
  if (mJutziStateUpdate) {
    state.enroll(*mfDSptr, make_policy<ReplaceBoundedState<Dimension, Scalar>>(1.0));
  } else {
    optionalOverridePolicy(SolidFieldNames::yieldStrength, make_policy<PorousYieldStrengthPolicy<Dimension>>());
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorosityModel<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDalphaDt);

  // Optional Jutzi deviatoric stress modifier
  if (mJutziStateUpdate) {
    derivs.enroll(*mfDSnewPtr);
  }
}

//------------------------------------------------------------------------------
// One time initializations at problem set up.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorosityModel<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Initialize the distention field.
  mAlpha = mAlpha0;

  // Get some state from the DataBase we're gonna need
  const auto  rhoFL = dataBase.fluidMassDensity();
  const auto& rho = **rhoFL.fieldForNodeList(mNodeList);

  // Solid density
  mSolidMassDensity = mAlpha0*rho;
  mSolidMassDensity.name(SolidFieldNames::porositySolidDensity);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorosityModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mMaxAbsDalphaDt, pathName + "/maxAbsDalphaDt");
  file.write(mAlpha0, pathName + "/alpha0");
  file.write(mAlpha, pathName + "/alpha");
  file.write(mDalphaDt, pathName + "/DalphaDt");
  file.write(mSolidMassDensity, pathName + "/solidMassDensity");
  if (mJutziStateUpdate) {
    file.write(*mfDSptr, pathName + "/fDS");
    file.write(*mfDSnewPtr, pathName + "/fDSnew");
  }
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorosityModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mMaxAbsDalphaDt, pathName + "/maxAbsDalphaDt");
  file.read(mAlpha0, pathName + "/alpha0");
  file.read(mAlpha, pathName + "/alpha");
  file.read(mDalphaDt, pathName + "/DalphaDt");
  file.read(mSolidMassDensity, pathName + "/solidMassDensity");
  if (mJutziStateUpdate) {
    file.read(*mfDSptr, pathName + "/fDS");
    file.read(*mfDSnewPtr, pathName + "/fDSnew");
  }
}

//------------------------------------------------------------------------------
// Compute the current porosity
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
PorosityModel<Dimension>::
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
