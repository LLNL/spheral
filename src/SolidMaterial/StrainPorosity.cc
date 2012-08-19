//---------------------------------Spheral++----------------------------------//
// StrainPorosity -- Strain porosity EOS modifier.
// 
// See header for references and such.
//----------------------------------------------------------------------------//

#include "StrainPorosity.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;
using NodeSpace::NodeList;
using PhysicsSpace::Physics;
using DataBaseSpace::DataBase;
using FileIOSpace::FileIO;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrainPorosity<Dimension>::
StrainPorosity(PorousEquationOfState<Dimension>& porousEOS,
               const NodeList<Dimension>& nodeList,
               const double phi0,
               const double epsE,
               const double epsX,
               const double kappa):
  Physics<Dimension>(),
  mAlpha0(1.0/(1.0 - phi0)),
  mEpsE(epsE),
  mEpsX(epsX),
  mKappa(kappa),
  mEpsC(0.0),
  mPorousEOS(porousEOS),
  mNodeList(nodeList),
  mAlpha(SolidFieldNames::porosityAlpha, nodeList, 1.0/(1.0 - phi0)),
  mDalphaDt(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, nodeList),
  mStrain(SolidFieldNames::porosityStrain, nodeList),
  mDstrainDt(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, nodeList),
  mRestart(DataOutput::registerWithRestart(*this)) {
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

  // TAU timers.
  TAU_PROFILE("StrainPorosity", "::evaluateDerivatives", TAU_USER);

  // Get the state fields.
  typedef typename State<Dimension>::KeyType Key;
  const Key gradvKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, mNodeList.name());
  const Key strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityStrain, mNodeList.name());
  const Key DalphaDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, mNodeList.name());
  const Key DstrainDtKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, mNodeList.name());
  const Field<Dimension, Tensor>& DvDx = derivs.field(gradvKey, Tensor::zero);
  const Field<Dimension, Scalar>& strain = state.field(strainKey, 0.0);
  Field<Dimension, Scalar>& DstrainDt = derivs.field(DstrainDtKey, 0.0);
  Field<Dimension, Scalar>& DalphaDt = derivs.field(DalphaDtKey, 0.0);

  const double A = mKappa*mAlpha0;
  const Scalar B = 2.0*(1.0 - mAlpha0*exp(mKappa*(mEpsX - mEpsE)))/FastMath::square(mEpsC - mEpsX);

  // Walk the nodes.
  for (unsigned i = 0; i != mNodeList.numInternalNodes(); ++i) {
    DstrainDt(i) = DvDx(i).Trace();
    const Scalar strainRate = min(0.0, DstrainDt(i));

    // The time rate of change for alpha depends on if we're in the exponential or
    // power-law regime of compaction.
    if (strain(i) <= mEpsE) {
      if (strain(i) >= mEpsX) {
        // Exponential.
        DalphaDt(i) = min(0.0, A*exp(mKappa*(strain(i) - mEpsE))*strainRate);
      } else {
        // Power-law.
        DalphaDt(i) = min(0.0, B*(mEpsE - strain(i))*strainRate);
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

  // TAU timers.
  TAU_PROFILE("StrainPorosity", "::registerState", TAU_USER);

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

  // TAU timers.
  TAU_PROFILE("StrainPorosity", "::registerDerivatives", TAU_USER);

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
}

