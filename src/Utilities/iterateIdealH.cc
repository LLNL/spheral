//------------------------------------------------------------------------------
// Iterate the ideal H algorithm to converge on a new H field.
// This routine replaces the H field in place.
//------------------------------------------------------------------------------
#include "iterateIdealH.hh"
#include "Field/FieldList.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/range.hh"
#include "Distributed/Communicator.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Geometry/GeometryRegistrar.hh"

#include <ctime>
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

template<typename Dimension>
void
iterateIdealH(DataBase<Dimension>& dataBase,
              SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              const vector<Boundary<Dimension>*>& boundaries,
              const int maxIterations,
              const double tolerance,
              const double nPerhForIteration,
              const bool sphericalStart,
              const bool fixDeterminant) {

  using SymTensor = typename Dimension::SymTensor;

  // Start the timing.
  const auto t0 = clock();

  // Extract the state we care about.
  const auto pos = dataBase.fluidPosition();
  auto m = dataBase.fluidMass();
  auto rho = dataBase.fluidMassDensity();
  auto H = dataBase.fluidHfield();

  // If we're using the fixDeterminant, take a snapshot of the input H determinants.
  auto Hdet0 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet0(itr) = H(itr).Determinant();
  }

  // Store the input nperh for each NodeList.
  // If we're rescaling the nodes per h for our work, make a cut at it.
  vector<double> nperh0;
  // Pulled divide by nPerhForIteration out of loop to improve optimization
  if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      const auto nperh = nodeListPtr->nodesPerSmoothingScale();
      nperh0.push_back(nperh);
      auto& Hfield = **(H.fieldForNodeList(*nodeListPtr));
      Hfield *= Dimension::rootnu(nperh / nPerhForIteration);
      nodeListPtr->nodesPerSmoothingScale(nPerhForIteration);
    }
  }
  else {
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      const auto nperh = nodeListPtr->nodesPerSmoothingScale();
      nperh0.push_back(nperh);
    }
  }
  CHECK(nperh0.size() == dataBase.numFluidNodeLists());

  // If we are both fixing the H determinant and rescaling the nodes per h,
  // we need a snapshot of the rescaled H determinant as well.
  auto Hdet1 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet1(itr) = H(itr).Determinant();
  }

  // If requested, start by making all the H's round (volume preserving).
  if (sphericalStart) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      const auto Hdeti = H(itr).Determinant();
      const auto Hi = Dimension::rootnu(Hdeti) * SymTensor::one;
      H(itr) = Hi;
    }
  }

  // Build a list of flags to indicate which nodes have been completed.
  auto flagNodeDone = dataBase.newFluidFieldList(0, "node completed");

  // Prepare the state and derivatives
  vector<Physics<Dimension>*> packages = {&smoothingScaleMethod};
  State<Dimension> state(dataBase, packages);
  StateDerivatives<Dimension> derivs(dataBase, packages);

  // Iterate until we either hit the max iterations or the H's achieve convergence.
  const auto numNodeLists = dataBase.numFluidNodeLists();
  auto maxDeltaH = 2.0*tolerance;
  auto itr = 0;
  while (itr < maxIterations and maxDeltaH > tolerance) {
    ++itr;
    maxDeltaH = 0.0;
    // flagNodeDone = 0;

    // Remove any old ghost node information from the NodeLists.
    for (auto k = 0u; k < numNodeLists; ++k) {
      auto nodeListPtr = *(dataBase.fluidNodeListBegin() + k);
      nodeListPtr->numGhostNodes(0);
      nodeListPtr->neighbor().updateNodes();
    }

    // Enforce boundary conditions.
    for (auto k = 0u; k < boundaries.size(); ++k) {
      auto boundaryPtr = *(boundaries.begin() + k);
      boundaryPtr->setAllGhostNodes(dataBase);
      boundaryPtr->applyFieldListGhostBoundary(m);
      boundaryPtr->applyFieldListGhostBoundary(rho);
      boundaryPtr->finalizeGhostBoundary();
      for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) nodeListPtr->neighbor().updateNodes();
    }

    // Call the smoothing scale package to get a new vote on the ideal H
    smoothingScaleMethod.initialize(0.0, 1.0, dataBase, state, derivs);
    derivs.Zero();
    smoothingScaleMethod.evaluateDerivatives(0.0, 1.0, dataBase, state, derivs);
    smoothingScaleMethod.finalizeDerivatives(0.0, 1.0, dataBase, state, derivs);
    
    // Extract the new ideal H vote
    auto H1 = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);

    // Set the new H and measure how much it changed
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto nodeListPtr = *(dataBase.fluidNodeListBegin() + nodeListi);
      const auto ni = nodeListPtr->numInternalNodes();

#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        if (flagNodeDone(nodeListi, i) == 0) {

          // If we are preserving the determinant, do it.
          if (fixDeterminant) {
            H1(nodeListi, i) *= Dimension::rootnu(Hdet1(nodeListi, i)/H1(nodeListi, i).Determinant());
            CHECK(fuzzyEqual(H1(nodeListi, i).Determinant(), Hdet1(nodeListi, i)));
          }

          // Check how much this H has changed.
          const auto H1sqrt = H1(nodeListi, i).sqrt();
          const auto phi = (H1sqrt*H(nodeListi, i).Inverse()*H1sqrt).Symmetric().eigenValues();
          const auto phimin = phi.minElement();
          const auto phimax = phi.maxElement();
          const auto deltaHi = max(abs(phimin - 1.0), abs(phimax - 1.0));
          if (deltaHi <= tolerance) flagNodeDone(nodeListi, i) = 1;
          maxDeltaH = max(maxDeltaH, deltaHi);

          // Assign the new H
          H(nodeListi, i) = H1(nodeListi, i);
        }
      }
    }

    // Globally reduce the max H change.
    maxDeltaH = allReduce(maxDeltaH, MPI_MAX, Communicator::communicator());

    // Output the statitics.
    if (Process::getRank() == 0)
      cerr << "iterateIdealH: (iteration, deltaH) = ("
           << itr << ", "
           << maxDeltaH << ")"
           << endl;

  }

  // If we have rescaled the nodes per h, now we have to iterate the H determinant
  // to convergence with the proper nperh.
  if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {

    // Reset the nperh.
    for (auto [k, nodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      CHECK(k < nperh0.size());
      //const double nperh = nperh0[k];
//       Field<Dimension, SymTensor>& Hfield = **(H.fieldForNodeList(**nodeListItr));
//       Hfield *= Dimension::rootnu(nPerhForIteration/nperh);
      nodeListPtr->nodesPerSmoothingScale(nperh0[k]);
    }
  }

  // If we're fixing the determinant, restore them.
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      H(itr) *= Dimension::rootnu(Hdet0(itr)/H(itr).Determinant());
      ENSURE(fuzzyEqual(H(itr).Determinant(), Hdet0(itr), 1.0e-10));
    }
  }

  // Leave the boundary conditions properly enforced.
  for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    nodeListPtr->numGhostNodes(0);
    nodeListPtr->neighbor().updateNodes();
  }
  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) {
    boundaryPtr->setAllGhostNodes(dataBase);
    boundaryPtr->finalizeGhostBoundary();
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      nodeListPtr->neighbor().updateNodes();
    }
  }

  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) boundaryPtr->applyFieldListGhostBoundary(m);
  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) boundaryPtr->finalizeGhostBoundary();

  // Report the final timing.
  const auto t1 = clock();
  if (Process::getRank() == 0)
    cerr << "iterateIdealH: required a total of "
         << (t1 - t0)/CLOCKS_PER_SEC
         << " seconds."
         << endl;
}

}

