//------------------------------------------------------------------------------
// Iterate the ideal H algorithm to converge on a new H field.
// This routine replaces the H field in place.
//------------------------------------------------------------------------------
#include <ctime>
#include "iterateIdealH.hh"
#include "Field/FieldList.hh"
#include "NodeList/SmoothingScaleBase.hh"

#ifdef USE_MPI
#include "mpi.h"
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {

using namespace std;

using DataBaseSpace::DataBase;
using BoundarySpace::Boundary;
using KernelSpace::TableKernel;
using FieldSpace::FieldList;
using FieldSpace::Field;
using NeighborSpace::ConnectivityMap;
using NodeSpace::SmoothingScaleBase;

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
  const unsigned rank = Process::getRank();

  // Start the timing.
  const clock_t t0 = clock();

  // Extract the state we care about.
  const FieldList<Dimension, Vector> r = dataBase.fluidPosition();
  FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();

  // If we're using the fixDeterminant, take a snapshot of the input H determinants.
  FieldList<Dimension, double> Hdet0 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (InternalNodeIterator<Dimension> itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet0(itr) = H(itr).Determinant();
  }

  // Store the input nperh for each NodeList.
  // If we're rescaling the nodes per h for our work, make a cut at it.
  vector<double> nperh0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd(); 
       ++nodeListItr) {
    const double nperh = (*nodeListItr)->nodesPerSmoothingScale();
    nperh0.push_back(nperh);
    if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {
      Field<Dimension, SymTensor>& Hfield = **(H.fieldForNodeList(**nodeListItr));
      Hfield *= Dimension::rootnu(nperh/nPerhForIteration);
      (*nodeListItr)->nodesPerSmoothingScale(nPerhForIteration);
    }
  }
  CHECK(nperh0.size() == dataBase.numFluidNodeLists());

  // If we are both fixing the H determinant and rescaling the nodes per h,
  // we need a snapshot of the rescaled H determinant as well.
  FieldList<Dimension, double> Hdet1 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (InternalNodeIterator<Dimension> itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet1(itr) = H(itr).Determinant();
  }

  // If requested, start by making all the H's round (volume preserving).
  if (sphericalStart) {
    for (InternalNodeIterator<Dimension> itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      const double Hdeti = H(itr).Determinant();
      const SymTensor Hi = Dimension::rootnu(Hdeti) * SymTensor::one;
      H(itr) = Hi;
    }
  }

  // Keep track of the step-wise changes in the H.
  FieldList<Dimension, double> deltaH = dataBase.newFluidFieldList(double());
  deltaH = 2.0*tolerance;

  // Iterate until we either hit the max iterations or the H's achieve convergence.
  const int numNodeLists = dataBase.numFluidNodeLists();
  int itr = 0;
  double maxDeltaH = 2.0*tolerance;
  while (itr < maxIterations && maxDeltaH > tolerance) {
    ++itr;
    maxDeltaH = 0.0;

    // Remove any old ghost node information from the NodeLists.
    for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd(); 
         ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
      (*nodeListItr)->neighbor().updateNodes();
    }

    // Enforce boundary conditions.
    for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
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

    // Build a list of flags to indicate which nodes have been completed.
    FieldList<Dimension, int> flagNodeDone = dataBase.newFluidFieldList(int());
    flagNodeDone = 0;

    // Any nodes that have already converged we flag as done.
    for (InternalNodeIterator<Dimension> nodeItr = dataBase.fluidInternalNodeBegin();
         nodeItr != dataBase.fluidInternalNodeEnd();
         ++nodeItr) {
      if (deltaH(nodeItr) <= tolerance) flagNodeDone(nodeItr) = 1;
    }

    // Prepare a FieldList to hold the new H.
    FieldList<Dimension, SymTensor> H1(H);
    H1.copyFields();

    // Get the new connectivity.
    dataBase.updateConnectivityMap();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();

    // Iterate over the NodeLists.
    int nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd(); 
         ++nodeListItr, ++nodeListi) {
      const Scalar hmin = (**nodeListItr).hmin();
      const Scalar hmax = (**nodeListItr).hmax();
      const Scalar hminratio = (**nodeListItr).hminratio();
      const Scalar nPerh = (**nodeListItr).nodesPerSmoothingScale();
      const int maxNumNeighbors = (**nodeListItr).maxNumNeighbors();

      // Iterate over the internal nodes of this NodeList.
      for (int i = 0; i != (**nodeListItr).numInternalNodes(); ++i) {

        // Has this node been done yet?
        if (flagNodeDone(nodeListi, i) == 0) {

          // Get the state and neighbors for this node.
          const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(*nodeListItr, i);
          const Vector& ri = r(nodeListi, i);
          const SymTensor& Hi = H(nodeListi, i);

          // Prepare to accumulate the zeroth and second moments for this node.
          Scalar zerothMoment = 0.0;
          SymTensor secondMoment;
          int numNeighbors = 0;

          // Iterate over the neighbor NodeLists.
          for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
            const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);

            // Neighbors from this NodeList.
            const vector<int>& connectivity = fullConnectivity[nodeListj];
            for (vector<int>::const_iterator jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {
              const int j = *jItr;

              // Increment the moments.
              const Vector& rj = r(nodeListj, j);
              const Vector rij = ri - rj;
              const Scalar etai = (Hi*rij).magnitude();
              const Scalar Wi = std::abs(W.gradValue(etai, 1.0));
              const SymTensor thpt = rij.selfdyad()/(rij.magnitude2() + 1.0e-10);
              zerothMoment += fweightij*Wi;
              secondMoment += fweightij*FastMath::square(Wi/Dimension::pownu1(etai + 1.0e-10))*thpt;
              ++numNeighbors;

            }
          }

          // Finish the moments and measure the new H.
          zerothMoment = Dimension::rootnu(zerothMoment);
          H1(nodeListi, i) = smoothingScaleMethod.newSmoothingScale(Hi,
                                                                    ri,
                                                                    zerothMoment,
                                                                    secondMoment,
                                                                    numNeighbors,
                                                                    W,
                                                                    hmin,
                                                                    hmax,
                                                                    hminratio,
                                                                    nPerh,
                                                                    maxNumNeighbors);

          // If we are preserving the determinant, do it.
          if (fixDeterminant) {
            H1(nodeListi, i) *= Dimension::rootnu(Hdet1(nodeListi, i)/H1(nodeListi, i).Determinant());
            CHECK(fuzzyEqual(H1(nodeListi, i).Determinant(), Hdet1(nodeListi, i)));
          }

          // Check how much this H has changed.
          const SymTensor H1sqrt = H1(nodeListi, i).sqrt();
          const Vector phi = (H1sqrt*Hi.Inverse()*H1sqrt).Symmetric().eigenValues();
          const double phimin = phi.minElement();
          const double phimax = phi.maxElement();
          const double deltaHi = max(abs(phimin - 1.0), abs(phimax - 1.0));
          deltaH(nodeListi, i) = deltaHi;
          maxDeltaH = max(maxDeltaH, deltaHi);
        }

        // Flag this node as completed.
        flagNodeDone(nodeListi, i) = 1;
      }
    }

    BEGIN_CONTRACT_SCOPE;
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
    END_CONTRACT_SCOPE;

    // Assign the new H's.
    H.assignFields(H1);

#ifdef USE_MPI
    {
      // Globally reduce the max H change.
      double tmp = maxDeltaH;
      MPI_Allreduce(&tmp, &maxDeltaH, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
    }
#endif

    // Output the statitics.
#ifdef USE_MPI
    if (rank == 0)
#endif
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
    for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
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
    for (InternalNodeIterator<Dimension> itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      H(itr) *= Dimension::rootnu(Hdet0(itr)/H(itr).Determinant());
      ENSURE(fuzzyEqual(H(itr).Determinant(), Hdet0(itr), 1.0e-10));
    }
  }


  // Leave the boundary conditions properly enforced.
  for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd(); 
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }
  for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
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

  FieldList<Dimension, Scalar> m = dataBase.fluidMass();
  for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
       boundaryItr != boundaries.end();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(m);
  }
  for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
       boundaryItr != boundaries.end();
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }

  // Report the final timing.
  const clock_t t1 = clock();
#ifdef USE_MPI
    if (rank == 0)
#endif
    cerr << "iterateIdealH: required a total of "
         << (t1 - t0)/CLOCKS_PER_SEC
         << " seconds."
         << endl;
}

}

