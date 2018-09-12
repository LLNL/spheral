//---------------------------------Spheral++----------------------------------//
// StrainPorosity -- Strain porosity EOS modifier.
// 
// See header for references and such.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/IncrementFieldList.hh"
#include "StrainPorosity.hh"

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
StrainPorosity<Dimension>::
StrainPorosity(PorousEquationOfState<Dimension>& porousEOS,
               PorousStrengthModel<Dimension>& porousStrength,
               const NodeList<Dimension>& nodeList,
               const double phi0,
               const double epsE,
               const double epsX,
               const double kappa,
               const double gammaS0,
               const double cS0,
               const double c0):
  Physics<Dimension>(),
  mAlpha0(1.0/(1.0 - phi0)),
  mEpsE(epsE),
  mEpsX(epsX),
  mKappa(kappa),
  mEpsC(0.0),
  mGammaS0(gammaS0),
  mcS0(cS0),
  mc0(c0),
  mPorousEOS(porousEOS),
  mPorousStrength(porousStrength),
  mNodeList(nodeList),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList, 1.0/(1.0 - phi0)),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mStrain(SolidFieldNames::porosityStrain, nodeList),
  mDstrainDt(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, nodeList),
  mRestart(registerWithRestart(*this)) {
  VERIFY2(mEpsE <= 0.0,
          "ERROR : epsE required to be epsE <= 0.0.");
  VERIFY2(mEpsX <= mEpsE,
          "StrainPorosity ERROR : epsX required to be epsX <= epsE.");
  VERIFY2(phi0 >= 0.0 and phi0 < 1.0,
          "ERROR : Initial porosity required to be in the range phi0 = [0.0, 1.0)");
  VERIFY2(kappa >= 0.0 and kappa <= 1.0,
          "ERROR : kappa required to be in range kappa = [0.0, 1.0]");
  mEpsC = 2.0*(1.0 - mAlpha0*exp(mKappa*(mEpsX - mEpsE)))/(mKappa*mAlpha0*exp(mKappa*(mEpsX - mEpsE))) + mEpsX;
  mPorousEOS.alpha(mAlpha);
  mPorousEOS.alpha0(mAlpha0);
  mPorousEOS.c0(c0);
  mPorousStrength.alpha(mAlpha);
  ENSURE(mPorousEOS.valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrainPorosity<Dimension>::
~StrainPorosity() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the alpha field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state fields.
  typedef typename State<Dimension>::KeyType Key;
  const Key gradvKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, mNodeList.name());
  const Key strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityStrain, mNodeList.name());
  const Key alphaKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, mNodeList.name());
  const Key DalphaDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, mNodeList.name());
  const Key DstrainDtKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, mNodeList.name());
  const Key DuDtKey = State<Dimension>::buildFieldKey(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, mNodeList.name());
  const Field<Dimension, Tensor>& DvDx = derivs.field(gradvKey, Tensor::zero);
  const Field<Dimension, Scalar>& strain = state.field(strainKey, 0.0);
  const Field<Dimension, Scalar>& alpha = state.field(alphaKey, 0.0);
  const Field<Dimension, Scalar>& DuDt = derivs.field(DuDtKey, 0.0);
  Field<Dimension, Scalar>& DstrainDt = derivs.field(DstrainDtKey, 0.0);
  Field<Dimension, Scalar>& DalphaDt = derivs.field(DalphaDtKey, 0.0);

  const Scalar PLcoef = 2.0*(1.0 - mAlpha0*exp(mKappa*(mEpsX - mEpsE)))/FastMath::square(mEpsC - mEpsX);

  // Walk the nodes.
  for (unsigned i = 0; i != mNodeList.numInternalNodes(); ++i) {
    DstrainDt(i) = DvDx(i).Trace();
    const Scalar alphai = alpha(i);
    const Scalar epsi = strain(i);
    const Scalar DuDti = DuDt(i);
    const Scalar strainRate = min(0.0, DstrainDt(i) - mGammaS0*safeInv(mcS0*mcS0)*DuDti);

    // The time rate of change for alpha depends on if we're in the elastic, exponential, or
    // power-law regime of compaction.
    if (epsi < 0.0) {
      if (epsi >= mEpsE) {
        // Elastic regime
        // We assume compaction is reversible in this regime.
        const double cs = mcS0 + (alphai - 1.0)*safeInv(mAlpha0 - 1.0)*(mc0 - mcS0);
        DalphaDt(i) = alphai*(1.0 - FastMath::square(cs*safeInv(mcS0)))*strainRate;
      } else if (epsi >= mEpsX) {
        // Exponential regime
        DalphaDt(i) = min(0.0, mKappa*alphai*strainRate);
      } else {
        // Power-law -- irreversible compaction.
        DalphaDt(i) = min(0.0, PLcoef*(mEpsE - epsi)*strainRate);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename StrainPorosity<Dimension>::TimeStepType
StrainPorosity<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return TimeStepType(1.0e100, "Rate of porosity change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;
  PolicyPointer strainPolicy(new IncrementState<Dimension, Scalar>());
  PolicyPointer alphaPolicy(new IncrementBoundedState<Dimension, Scalar, Scalar>(1.0));
  state.enroll(mStrain, strainPolicy);
  state.enroll(mAlpha, alphaPolicy);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDstrainDt);
  derivs.enroll(mDalphaDt);
}

//------------------------------------------------------------------------------
// One time initializations at problem set up.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Initialize the distention field.
  mAlpha = mAlpha0;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mAlpha, pathName + "/alpha");
  file.write(mDalphaDt, pathName + "/DalphaDt");
  file.write(mStrain, pathName + "/strain");
  file.write(mDstrainDt, pathName + "/DstrainDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mAlpha, pathName + "/alpha");
  file.read(mDalphaDt, pathName + "/DalphaDt");
  file.read(mStrain, pathName + "/strain");
  file.read(mDstrainDt, pathName + "/DstrainDt");
}

}
