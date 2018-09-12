//------------------------------------------------------------------------------
// Iterate the ideal H algorithm to converge on a new H field.
// This routine replaces the H field in place.
//------------------------------------------------------------------------------
#include "iterateIdealH.hh"
#include "Field/FieldList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Utilities/allReduce.hh"
#include "Distributed/Communicator.hh"

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
              const vector<Boundary<Dimension>*>& boundaries,
              const TableKernel<Dimension>& W,
              const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              const int maxIterations,
              const double tolerance,
              const double nPerhForIteration,
              const bool sphericalStart,
              const bool fixDeterminant) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Get the local rank.
  const auto rank = Process::getRank();

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
  for (auto nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd(); 
       ++nodeListItr) {
    const auto nperh = (*nodeListItr)->nodesPerSmoothingScale();
    nperh0.push_back(nperh);
    if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {
      auto& Hfield = **(H.fieldForNodeList(**nodeListItr));
      Hfield *= Dimension::rootnu(nperh/nPerhForIteration);
      (*nodeListItr)->nodesPerSmoothingScale(nPerhForIteration);
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

  // Keep track of the step-wise changes in the H.
  auto deltaH = dataBase.newFluidFieldList(double());
  deltaH = 2.0*tolerance;

  // Iterate until we either hit the max iterations or the H's achieve convergence.
  const auto numNodeLists = dataBase.numFluidNodeLists();
  auto maxDeltaH = 2.0*tolerance;
  int itr = 0;
  while (itr < maxIterations && maxDeltaH > tolerance) {
    ++itr;
    maxDeltaH = 0.0;

    // Remove any old ghost node information from the NodeLists.
    for (auto k = 0; k < numNodeLists; ++k) {
      auto nodeListPtr = *(dataBase.fluidNodeListBegin() + k);
      nodeListPtr->numGhostNodes(0);
      nodeListPtr->neighbor().updateNodes();
    }

    // Enforce boundary conditions.
    for (auto k = 0; k < boundaries.size(); ++k) {
      auto boundaryPtr = *(boundaries.begin() + k);
      boundaryPtr->setAllGhostNodes(dataBase);
      boundaryPtr->applyFieldListGhostBoundary(m);
      boundaryPtr->applyFieldListGhostBoundary(rho);
      boundaryPtr->finalizeGhostBoundary();
      for (auto nodeListItr = dataBase.fluidNodeListBegin();
           nodeListItr != dataBase.fluidNodeListEnd(); 
           ++nodeListItr) {
        (*nodeListItr)->neighbor().updateNodes();
      }
    }

    // Build a list of flags to indicate which nodes have been completed.
    auto flagNodeDone = dataBase.newFluidFieldList(int());
    flagNodeDone = 0;

    // Any nodes that have already converged we flag as done.
    for (auto nodeItr = dataBase.fluidInternalNodeBegin();
         nodeItr != dataBase.fluidInternalNodeEnd();
         ++nodeItr) {
      if (deltaH(nodeItr) <= tolerance) flagNodeDone(nodeItr) = 1;
    }

    // Prepare a FieldList to hold the new H.
    FieldList<Dimension, SymTensor> H1(H);
    H1.copyFields();

    // Get the new connectivity.
    dataBase.updateConnectivityMap(false);
    const auto& connectivityMap = dataBase.connectivityMap();

    // Iterate over the NodeLists.
    int nodeListi = 0;
    for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
      const auto nodeListPtr = *(dataBase.fluidNodeListBegin() + nodeListi);
      const Scalar hmin = nodeListPtr->hmin();
      const Scalar hmax = nodeListPtr->hmax();
      const Scalar hminratio = nodeListPtr->hminratio();
      const Scalar nPerh = nodeListPtr->nodesPerSmoothingScale();

      // Iterate over the internal nodes of this NodeList.
#pragma omp parallel
      {
        auto maxDeltaH_local = maxDeltaH;
#pragma omp for
        for (auto i = 0; i < nodeListPtr->numInternalNodes(); ++i) {

          // Has this node been done yet?
          if (flagNodeDone(nodeListi, i) == 0) {

            // Get the state and neighbors for this node.
            const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListPtr, i);
            const auto& posi = pos(nodeListi, i);
            const auto& Hi = H(nodeListi, i);
            const auto  mi = m(nodeListi, i);
            const auto  rhoi = rho(nodeListi, i);

            // Prepare to accumulate the zeroth and second moments for this node.
            auto zerothMoment = 0.0;
            SymTensor secondMoment;
            Scalar fweightij;

            // Iterate over the neighbor NodeLists.
            for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

              // Neighbors from this NodeList.
              const auto& connectivity = fullConnectivity[nodeListj];
              for (auto jItr = connectivity.begin();
                   jItr != connectivity.end();
                   ++jItr) {
                const int j = *jItr;

                // Increment the moments.
                const auto& posj = pos(nodeListj, j);
                const auto  mj = m(nodeListj, j);
                const auto  rhoj = rho(nodeListj, j);

                fweightij = 1.0;
                if (nodeListi != nodeListj) {
                  if (dataBase.isRZ) {
                    const auto ri = abs(posi.y());
                    const auto rj = abs(posj.y());
                    const auto mRZi = mi/(2.0*M_PI*ri);
                    const auto mRZj = mj/(2.0*M_PI*rj);
                    fweightij = mRZj*rhoi/(mRZi*rhoj);
                  } else {
                    fweightij = mj*rhoi/(mi*rhoj);
                  }
                }
                                        
                const auto xij = posi - posj;
                const auto etai = (Hi*xij).magnitude();
                const auto Wi = std::abs(W.gradValue(etai, 1.0));
                const auto thpt = xij.selfdyad()/(xij.magnitude2() + 1.0e-10);
                zerothMoment += fweightij*Wi;
                secondMoment += fweightij*FastMath::square(Wi*safeInvVar(xij.magnitude2()))*thpt;
                // secondMoment += fweightij*FastMath::square(Wi/Dimension::pownu1(etai + 1.0e-10))*thpt;
              }
            }

            // Finish the moments and measure the new H.
            zerothMoment = Dimension::rootnu(zerothMoment);
            H1(nodeListi, i) = smoothingScaleMethod.newSmoothingScale(Hi,
                                                                      posi,
                                                                      zerothMoment,
                                                                      secondMoment,
                                                                      W,
                                                                      hmin,
                                                                      hmax,
                                                                      hminratio,
                                                                      nPerh,
                                                                      connectivityMap,
                                                                      nodeListi,
                                                                      i);

            // If we are preserving the determinant, do it.
            if (fixDeterminant) {
              H1(nodeListi, i) *= Dimension::rootnu(Hdet1(nodeListi, i)/H1(nodeListi, i).Determinant());
              CHECK(fuzzyEqual(H1(nodeListi, i).Determinant(), Hdet1(nodeListi, i)));
            }

            // Check how much this H has changed.
            const auto H1sqrt = H1(nodeListi, i).sqrt();
            const auto phi = (H1sqrt*Hi.Inverse()*H1sqrt).Symmetric().eigenValues();
            const auto phimin = phi.minElement();
            const auto phimax = phi.maxElement();
            const auto deltaHi = max(abs(phimin - 1.0), abs(phimax - 1.0));
            deltaH(nodeListi, i) = deltaHi;
            maxDeltaH = max(maxDeltaH, deltaHi);
          }

          // Flag this node as completed.
          flagNodeDone(nodeListi, i) = 1;
        }

#pragma omp critical
        maxDeltaH = max(maxDeltaH, maxDeltaH_local);
      }
    }

    BEGIN_CONTRACT_SCOPE
    {
      // Ensure that all nodes have been calculated.
      for (typename FieldList<Dimension, int>::const_iterator fieldItr = flagNodeDone.begin();
           fieldItr != flagNodeDone.end();
           ++fieldItr) {
        for (int i = 0; i != (*fieldItr)->nodeListPtr()->numInternalNodes(); ++i) {
          CHECK((**fieldItr)[i] == 1);
        }
      }
    }
    END_CONTRACT_SCOPE

    // Assign the new H's.
    H.assignFields(H1);

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
    size_t k = 0;
    for (auto nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd(); 
         ++nodeListItr, ++k) {
      CHECK(k < nperh0.size());
      const double nperh = nperh0[k];
//       Field<Dimension, SymTensor>& Hfield = **(H.fieldForNodeList(**nodeListItr));
//       Hfield *= Dimension::rootnu(nPerhForIteration/nperh);
      (*nodeListItr)->nodesPerSmoothingScale(nperh0[k]);
    }
    CHECK(k == nperh0.size());

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
  for (auto nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd(); 
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }
  for (auto boundaryItr = boundaries.begin(); 
       boundaryItr != boundaries.end();
       ++boundaryItr) {
    (*boundaryItr)->setAllGhostNodes(dataBase);
    (*boundaryItr)->finalizeGhostBoundary();
    for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd(); 
         ++nodeListItr) {
      (*nodeListItr)->neighbor().updateNodes();
    }
  }

  for (auto boundaryItr = boundaries.begin(); 
       boundaryItr != boundaries.end();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(m);
  }
  for (auto boundaryItr = boundaries.begin(); 
       boundaryItr != boundaries.end();
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }

  // Report the final timing.
  const auto t1 = clock();
  if (Process::getRank() == 0)
    cerr << "iterateIdealH: required a total of "
         << (t1 - t0)/CLOCKS_PER_SEC
         << " seconds."
         << endl;
}

}

