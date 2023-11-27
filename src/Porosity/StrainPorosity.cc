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
#include "Porosity/StrainPorosity.hh"

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
StrainPorosity(const SolidNodeList<Dimension>& nodeList,
               const double phi0,
               const double epsE,
               const double epsX,
               const double kappa,
               const double gammaS0,
               const double cS0,
               const double c0,
               const double rhoS0):
  PorosityModel<Dimension>(nodeList, phi0, cS0, c0, rhoS0, false),
  mEpsE(epsE),
  mEpsX(epsX),
  mKappa(kappa),
  mGammaS0(gammaS0),
  mStrain(SolidFieldNames::porosityStrain, nodeList),
  mDstrainDt(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, nodeList) {
  VERIFY2(mEpsE <= 0.0,
          "ERROR : epsE required to be epsE <= 0.0.");
  VERIFY2(mEpsX <= mEpsE,
          "StrainPorosity ERROR : epsX required to be epsX <= epsE.");
  VERIFY2(kappa >= 0.0 and kappa <= 1.0,
          "ERROR : kappa required to be in range kappa = [0.0, 1.0]");
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrainPorosity<Dimension>::
StrainPorosity(const SolidNodeList<Dimension>& nodeList,
               const Field<Dimension, Scalar>& phi0,
               const double epsE,
               const double epsX,
               const double kappa,
               const double gammaS0,
               const double cS0,
               const Field<Dimension, Scalar>& c0,
               const double rhoS0):
  PorosityModel<Dimension>(nodeList, phi0, cS0, c0, rhoS0, false),
  mEpsE(epsE),
  mEpsX(epsX),
  mKappa(kappa),
  mGammaS0(gammaS0),
  mStrain(SolidFieldNames::porosityStrain, nodeList),
  mDstrainDt(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, nodeList) {
  VERIFY2(mEpsE <= 0.0,
          "ERROR : epsE required to be epsE <= 0.0.");
  VERIFY2(mEpsX <= mEpsE,
          "StrainPorosity ERROR : epsX required to be epsX <= epsE.");
  VERIFY2(kappa >= 0.0 and kappa <= 1.0,
          "ERROR : kappa required to be in range kappa = [0.0, 1.0]");
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
evaluateDerivatives(const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state fields.
  const auto  gradvKey = State<Dimension>::buildFieldKey(HydroFieldNames::velocityGradient, mNodeList.name());
  const auto  strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityStrain, mNodeList.name());
  const auto  alphaKey = State<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, mNodeList.name());
  const auto  DalphaDtKey = State<Dimension>::buildFieldKey(IncrementBoundedState<Dimension, Scalar, Scalar>::prefix() + SolidFieldNames::porosityAlpha, mNodeList.name());
  const auto  DstrainDtKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityStrain, mNodeList.name());
  const auto  DuDtKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, mNodeList.name());
  const auto& DvDx = derivs.field(gradvKey, Tensor::zero);
  const auto& strain = state.field(strainKey, 0.0);
  const auto& alpha = state.field(alphaKey, 0.0);
  const auto& DuDt = derivs.field(DuDtKey, 0.0);
  auto&       DstrainDt = derivs.field(DstrainDtKey, 0.0);
  auto&       DalphaDt = derivs.field(DalphaDtKey, 0.0);

  // Walk the nodes.
  const auto n = mNodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    DstrainDt(i) = DvDx(i).Trace();
    const auto alphai = alpha(i);
    const auto epsi = strain(i);
    const auto DuDti = DuDt(i);
    const auto strainRate = min(0.0, DstrainDt(i) - mGammaS0*safeInv(mcS0*mcS0)*DuDti);
    const auto alpha0i = mAlpha0(i);
    const auto c0i = mc0(i);

    // The time rate of change for alpha depends on if we're in the elastic, exponential, or
    // power-law regime of compaction.
    if (epsi < 0.0) {
      if (epsi >= mEpsE) {
        // Elastic regime
        // We assume compaction is reversible in this regime.
        const auto cs = mcS0 + (alphai - 1.0)*safeInv(alpha0i - 1.0)*(c0i - mcS0);
        DalphaDt(i) = alphai*(1.0 - FastMath::square(cs*safeInv(mcS0)))*strainRate;
      } else if (epsi >= mEpsX) {
        // Exponential regime
        DalphaDt(i) = min(0.0, mKappa*alphai*strainRate);
      } else {
        // Power-law -- irreversible compaction.
        const auto epsCi = 2.0*(1.0 - alpha0i*exp(mKappa*(mEpsX - mEpsE)))/(mKappa*alpha0i*exp(mKappa*(mEpsX - mEpsE))) + mEpsX;
        const auto PLcoefi = 2.0*(1.0 - alpha0i*exp(mKappa*(mEpsX - mEpsE)))/FastMath::square(epsCi - mEpsX);
        DalphaDt(i) = min(0.0, PLcoefi*(mEpsE - epsi)*strainRate);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  PorosityModel<Dimension>::registerState(dataBase, state);
  state.enroll(mStrain, make_policy<IncrementState<Dimension, Scalar>>());
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  PorosityModel<Dimension>::registerDerivatives(dataBase, derivs);
  derivs.enroll(mDstrainDt);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPorosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  PorosityModel<Dimension>::dumpState(file, pathName);
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
  PorosityModel<Dimension>::restoreState(file, pathName);
  file.read(mStrain, pathName + "/strain");
  file.read(mDstrainDt, pathName + "/DstrainDt");
}

}
