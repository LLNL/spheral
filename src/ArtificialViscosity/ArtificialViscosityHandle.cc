//---------------------------------Spheral++----------------------------------//
// ArtificialViscosityHandle -- A base class for ArtficialViscosity that strips
// off the QPiType template parameter.  This makes a convenient handle to break
// that template parameter from spreading into classes that need to consume an
// ArtificialViscosity.
//
// Created by JMO, Fri Dec 13 10:06:12 PST 2024
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"

#include "ArtificialViscosityHandle.hh"

#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosityHandle<Dimension>::
ArtificialViscosityHandle(const Scalar Clinear,
                          const Scalar Cquadratic,
                          const TableKernel<Dimension>& kernel):
  Physics<Dimension>(),
  mClinear(Clinear),
  mCquadratic(Cquadratic),
  mBalsaraShearCorrection(false),
  mEpsilon2(1.0e-2),
  mNegligibleSoundSpeed(1e-10),
  mMaxViscousPressure(FieldStorageType::CopyFields),
  mRigorousVelocityGradient(false),
  mWT(kernel),
  mDvDx(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Vote on a timestep.  We assume the hydro accounts for our contribution to the
// timestep, so by default no-vote.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Physics<Dimension>::TimeStepType
ArtificialViscosityHandle<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar currentTime) const {
  return TimeStepType(1e100, "ArtificialViscosity taken care of by hydro -- no vote");
}

//------------------------------------------------------------------------------
// Register state
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  if (this->requireVelocityGradient() or
      this->balsaraShearCorrection()) state.enroll(mDvDx);
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mMaxViscousPressure);
}

//------------------------------------------------------------------------------
// Apply ghost boundaries
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  if (this->requireVelocityGradient() or
      this->balsaraShearCorrection()) {
    auto DvDx = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);
    for (auto* bcPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
      bcPtr->applyFieldListGhostBoundary(DvDx);
    }
  }
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  // Size our internal state fields
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);

  // Prepare the initial velocity gradient
  if (this->requireVelocityGradient() or
      this->balsaraShearCorrection()) {
    dataBase.resizeFluidFieldList(mM, Tensor::zero, "AV M correction field", false);
    dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::ArtificialViscosityVelocityGradient, false);
    this->updateVelocityGradient(dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Post-state update.  Update the velocity gradient if needed.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ArtificialViscosityHandle<Dimension>::
postStateUpdate(const Scalar t,
                const Scalar dt,
                const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
  if (this->requireVelocityGradient() or
      this->balsaraShearCorrection()) {
    dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::ArtificialViscosityVelocityGradient, false);
    this->updateVelocityGradient(dataBase, state, derivs);
    return true;
  } else {
    return false;
  }
}

//------------------------------------------------------------------------------
// Update our velocity gradient
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
updateVelocityGradient(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivs) {

  auto DvDx_AV = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);

  if (mRigorousVelocityGradient) {

    // Measure a linearly corrected velocity gradient based on the corrected SPH form.
    // Same thing used in our SPH hydro measurement of DvDx.
    const auto& W = this->kernel();
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto& pairs = connectivityMap.nodePairList();
    const auto  numNodeLists = dataBase.numFluidNodeLists();
    const auto  npairs = pairs.size();

    // Grab the state we need
    const auto mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);

    // Clear our current estimate
    DvDx_AV.Zero();
    mM.Zero();

#pragma omp parallel
    {
      // Thread private scratch variables
      int i, j, nodeListi, nodeListj;
      Scalar Hdeti, Hdetj, etaMagi, etaMagj;
      Vector rij, vij, etai, etaj, etaiUnit, etajUnit, gradWi, gradWj;
      Tensor QPiij, QPiji;

      typename SpheralThreads<Dimension>::FieldListStack threadStack;
      auto DvDx_thread = DvDx_AV.threadCopy(threadStack);
      auto M_thread = mM.threadCopy(threadStack);

#pragma omp for
      for (auto kk = 0u; kk < npairs; ++kk) {
        i = pairs[kk].i_node;
        j = pairs[kk].j_node;
        nodeListi = pairs[kk].i_list;
        nodeListj = pairs[kk].j_list;

        // Get the state for node i.
        const auto& ri = position(nodeListi, i);
        const auto& mi = mass(nodeListi, i);
        const auto& vi = velocity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        Hdeti = Hi.Determinant();
        CHECK(mi > 0.0);
        CHECK(Hdeti > 0.0);

        auto& DvDxi = DvDx_thread(nodeListi, i);
        auto& Mi = M_thread(nodeListi, i);

        // Get the state for node j
        const auto& rj = position(nodeListj, j);
        const auto& mj = mass(nodeListj, j);
        const auto& vj = velocity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        Hdetj = Hj.Determinant();
        CHECK(mj > 0.0);
        CHECK(Hdetj > 0.0);

        auto& DvDxj = DvDx_thread(nodeListj, j);
        auto& Mj = M_thread(nodeListj, j);

        // Node displacement.
        rij = ri - rj;
        etai = Hi*rij;
        etaj = Hj*rij;
        etaMagi = etai.magnitude();
        etaMagj = etaj.magnitude();
        etaiUnit = etai*safeInvVar(etaMagi);
        etajUnit = etaj*safeInvVar(etaMagj);
        CHECK(etaMagi >= 0.0);
        CHECK(etaMagj >= 0.0);

        // Symmetrized kernel weight and gradient.
        gradWi = W.grad(etaMagi, Hdeti) * Hi * etaiUnit;
        gradWj = W.grad(etaMagj, Hdetj) * Hj * etajUnit;

        // Velocity gradient.
        vij = vi - vj;
        DvDxi -= mj*vij.dyad(gradWi);
        DvDxj -= mi*vij.dyad(gradWj);

        // Linear gradient correction term.
        Mi -= mj*rij.dyad(gradWi);
        Mj -= mi*rij.dyad(gradWj);
      } // loop over pairs

      // Reduce the thread values to the master.
      threadReduceFieldLists<Dimension>(threadStack);

    }   // OpenMP parallel region

    // Finish up the gradient for each point.
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto ni = DvDx_AV[nodeListi]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {

        // Get the state for node i.
        const auto& mi = mass(nodeListi, i);
        const auto& rhoi = massDensity(nodeListi, i);
        const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
        CHECK(mi > 0.0);
        CHECK(rhoi > 0.0);

        auto& DvDxi = DvDx_AV(nodeListi, i);
        auto& Mi = mM(nodeListi, i);

        // Finish the gradient of the velocity.
        CHECK(rhoi > 0.0);
        if (std::abs(Mi.Determinant()) > 1.0e-10 and
            numNeighborsi > Dimension::pownu(2)) {
          Mi = Mi.Inverse();
          DvDxi = DvDxi*Mi;
        } else {
          DvDxi /= rhoi;
        }
      }
    }

  } else {

    // Just copy the hydro estimate for the velocity gradient. This option incurs being
    // a bit off in the time level of DvDx during an integration cycle but is cheaper.
    const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
    DvDx_AV.assignFields(DvDx);
    for (auto* fptr: DvDx_AV) fptr->name(HydroFieldNames::ArtificialViscosityVelocityGradient);

  }

}

//------------------------------------------------------------------------------
// Dump the current state of the Q to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  if (this->requireVelocityGradient() or this->balsaraShearCorrection()) {
    file.write(mDvDx, pathName + "/velocity_gradient");
  }
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityHandle<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  if (this->requireVelocityGradient() or this->balsaraShearCorrection()) {
    file.read(mDvDx, pathName + "/velocity_gradient");
  }
}

}
