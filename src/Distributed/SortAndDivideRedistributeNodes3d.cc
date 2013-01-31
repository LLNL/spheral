//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the SortAndDivideNeighbor algorithm.
//----------------------------------------------------------------------------//
#include <algorithm>
#include <vector>
#include <map>

#include <fstream>
#include <cstdlib>

#include "TAU.h"

#include "SortAndDivideRedistributeNodes3d.hh"
#include "DomainNode.hh"
#include "BoundingVolumeDistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/EigenStruct.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Communicator.hh"

#include "DBC.hh"
#include "cdebug.hh"

namespace Spheral {
namespace PartitionSpace {

using namespace std;

using DataBaseSpace::DataBase;
using NodeSpace::NodeList;
using BoundarySpace::BoundingVolumeDistributedBoundary;
using BoundarySpace::Boundary;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Construct with the given node extent.
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes3d::
SortAndDivideRedistributeNodes3d(const double Hextent):
  SortAndDivideRedistributeNodes<Dim<3> >(Hextent) {
  cdebug << "SortAndDivideRedistributeNodes3d::SortAndDivideRedistributeNodes3d()" << endl;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes3d::
~SortAndDivideRedistributeNodes3d() {
  cdebug << "SortAndDivideRedistributeNodes3d::~SortAndDivideRedistributeNodes3d()" << endl;
}

//------------------------------------------------------------------------------
// The main method of this class.  Repartition the node distribution.
//------------------------------------------------------------------------------
void
SortAndDivideRedistributeNodes3d::
redistributeNodes(DataBase<Dim<3> >& dataBase,
                  vector<Boundary<Dim<3> >*> boundaries) {
  cdebug << "SortAndDivideRedistributeNodes3d::redistributeNodes" << endl;

  // TAU timers.
  TAU_PROFILE("SortAndDivideRedistributeNodes3d", "::redistributeNodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3GlobalIDs, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : build global node IDs", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3CurrentDomain, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : build current domain decomposition", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3GhostNodes, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : build ghost nodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3LocalWork, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : measure local work", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3Shape, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : measure shape tensor", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3GlobalWork, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : measure global work", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3XSort, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Sort all nodes by x", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3YSort, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Sort slab nodes by y", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3ZSort, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Sort column nodes by z", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3XPop, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Pop x nodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3YPop, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Pop y nodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeSDR3ZPop, "SortAndDivideRedistributeNodes3d", "::redistributeNodes : Pop z nodes", TAU_USER);

  // Number of processors.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();

  // Go over each NodeList, and clear out any ghost nodes.
  cdebug << "Removing ghost nodes." << endl;
  for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // If the user did not specify any boundaries, then create a Distributed
  // boundary for local use.
  // Note that if boundary conditions were passed in, we assume that the Distributed
  // boundary is already in there.
  TAU_PROFILE_START(TimeSDR3GhostNodes);
  cdebug << "Building SortAndDivideDistributedBoundary" << endl;
  BoundingVolumeDistributedBoundary<Dimension>& bound = BoundingVolumeDistributedBoundary<Dimension>::instance();
  if (boundaries.size() == 0) boundaries.push_back(&bound);

  // Use the boundaries to create the ghost nodes.  We need this to compute the
  // work per node (defined as the number of neighbors) correctly.
  for (vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) {
    (*boundItr)->setAllGhostNodes(dataBase);
    (*boundItr)->finalizeGhostBoundary();
    for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
  }
  TAU_PROFILE_STOP(TimeSDR3GhostNodes);

  // Get the work per node.
  TAU_PROFILE_START(TimeSDR3LocalWork);
  const FieldList<Dimension, Scalar> work = this->workPerNode(dataBase, Hextent());
  TAU_PROFILE_STOP(TimeSDR3LocalWork);

  // Once again clear out any ghost nodes.
  for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Build the set of global node IDs.
  TAU_PROFILE_START(TimeSDR3GlobalIDs);
  const FieldList<Dimension, int> globalIDs = NodeSpace::globalNodeIDs(dataBase);
  TAU_PROFILE_STOP(TimeSDR3GlobalIDs);

  // Get the local description of the domain distribution.
  TAU_PROFILE_START(TimeSDR3CurrentDomain);
  vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs, work);
  TAU_PROFILE_STOP(TimeSDR3CurrentDomain);

  // Clear the domain assignments of all nodes.
  for (vector<DomainNode<Dimension> >::iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) itr->domainID = -1;

  // Compute the shape tensor describing the total node distribution in space.
  TAU_PROFILE_START(TimeSDR3Shape);
  const EigenStruct<3> shapeTensor = this->shapeTensor(nodeDistribution);
  TAU_PROFILE_STOP(TimeSDR3Shape);

  // Rotate the nodes into this frame.
  this->rotateIntoShapeTensorFrame(shapeTensor, nodeDistribution);

  // Get the number of domains we should have per constant work chunk in the primary direction.
  vector< vector<int> > domainsPerStep = domainsPerChunk(shapeTensor);

  // Output the initial load distribution statistics.
  const string initialLoadStats = this->gatherDomainDistributionStatistics(work);
  if (procID == 0) {
    cerr << "SortAndDivideRedistributeNodes::redistributeNodes initial load balancing:" << endl
         << initialLoadStats << endl
         << "    Domain distribution shape tensor: " << shapeTensor.eigenValues << endl;
    for (int i = 0; i != Dimension::nDim; ++i) {
      cerr << "    " << shapeTensor.eigenVectors.getColumn(i) << endl;
    }
  }
    
  // Compute the total work, and the target work per processor.
  TAU_PROFILE_START(TimeSDR3GlobalWork);
  double localWork = 0.0;
  for (vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) localWork += itr->work;
  CHECK(numProcs > 0);
  double globalWork;
  MPI_Allreduce(&localWork, &globalWork, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
  const double targetWorkPerDomain = globalWork / numProcs;
  CHECK(distinctlyGreaterThan(targetWorkPerDomain, 0.0));
  TAU_PROFILE_STOP(TimeSDR3GlobalWork);

  // Copy the domain nodes to a list.
  list<DomainNode<Dimension> > sortedNodeDistribution;
  for (vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) sortedNodeDistribution.push_back(*itr);

  // Sort the nodes by x position.
  TAU_PROFILE_START(TimeSDR3XSort);
  this->sortByPositions(sortedNodeDistribution, 0);
  TAU_PROFILE_STOP(TimeSDR3XSort);

  // The new node distribution we're going to build.
  vector<DomainNode<Dimension> > newNodeDistribution;
  newNodeDistribution.reserve(nodeDistribution.size());

  // Iterate over the x domain indicies.
  int assignDomainID = 0;
  for (int ix = 0; ix != domainsPerStep.size(); ++ix) {

    // Pop off the set of nodes for the y-z slab in this step.
    TAU_PROFILE_START(TimeSDR3XPop);
    const vector<int>& numYZChunks = domainsPerStep[ix];
    CHECK(numYZChunks.size() > 0);
    int numYZDomains = 0;
    for (vector<int>::const_iterator itr = numYZChunks.begin();
         itr != numYZChunks.end();
         ++itr) numYZDomains += *itr;
    const double targetYZChunkWork = numYZDomains*targetWorkPerDomain;

    list<DomainNode<Dimension> > yzchunkNodes;
    if (ix == domainsPerStep.size() - 1) {
      yzchunkNodes = sortedNodeDistribution;
    } else {
      yzchunkNodes = this->popFrontNodes(sortedNodeDistribution,
                                         targetYZChunkWork,
                                         0);
    }
    TAU_PROFILE_STOP(TimeSDR3XPop);

    // Re-sort this set of nodes by their y positions.
    TAU_PROFILE_START(TimeSDR3YSort);
    this->sortByPositions(yzchunkNodes, 1);
    TAU_PROFILE_STOP(TimeSDR3YSort);

    // Iterate over the number of y domains we're assigning for this slab of work.
    for (int iy = 0; iy != numYZChunks.size(); ++iy) {

      // Pop off a chunk of nodes we'll divvy up in the z direction.
      TAU_PROFILE_START(TimeSDR3YPop);
      const int numZChunks = numYZChunks[iy];
      CHECK(numZChunks > 0);
      const double targetZChunkWork = numZChunks*targetWorkPerDomain;

      // If this is the last domain in this slice, make sure we got all the nodes.
      // Otherwise, just pop off the target domains' work.
      list<DomainNode<Dimension> > chunkNodes;
      if (iy == numYZChunks.size() - 1) {
        chunkNodes = yzchunkNodes;
      } else {
        chunkNodes = this->popFrontNodes(yzchunkNodes,
                                         targetZChunkWork,
                                         1);
      }
      TAU_PROFILE_STOP(TimeSDR3YPop);

      // Re-sort this set of nodes by their z positions.
      TAU_PROFILE_START(TimeSDR3ZSort);
      this->sortByPositions(chunkNodes, 2);
      TAU_PROFILE_STOP(TimeSDR3ZSort);

      // Iterator over the number of z domains we'll be assigning.
      for (int iz = 0; iz != numZChunks; ++iz) {

        if (procID == 0) cerr << "Assigning domain " << assignDomainID 
                              << " of " << numProcs << "...";

        // Peel off nodes from the front of the unassigned nodes, until the desired work
        // load for this domain is reached.  Note that in this step we use the z index as
        // the primary sorting comparison.
        TAU_PROFILE_START(TimeSDR3ZPop);
        list<DomainNode<Dimension> > thisDomainNodes;
        if (iz == numZChunks - 1) {
          thisDomainNodes = chunkNodes;
        } else {
          thisDomainNodes = this->popFrontNodes(chunkNodes,
                                                targetWorkPerDomain,
                                                2);
        }
        TAU_PROFILE_STOP(TimeSDR3ZPop);

        // Assign the nodes to their new domain.
        CHECK(assignDomainID >= 0 && assignDomainID < numProcs);
        for (list<DomainNode<Dimension> >::const_iterator itr = thisDomainNodes.begin();
             itr != thisDomainNodes.end();
             ++itr) {
          newNodeDistribution.push_back(*itr);
          newNodeDistribution.back().domainID = assignDomainID;
        }

        // Increment the domain we're assigning to.
        ++assignDomainID;
        if (procID == 0) cerr << "Done." << endl;

      }

//       // Assign any remaining nodes from this chunks work to the last domain we were working
//       // on.
//       CHECK(assignDomainID - 1 >= 0 && assignDomainID - 1 < numProcs);
//       for (list<DomainNode<Dimension> >::const_iterator itr = chunkNodes.begin();
//            itr != chunkNodes.end();
//            ++itr) {
//         newNodeDistribution.push_back(*itr);
//         newNodeDistribution.back().domainID = assignDomainID - 1;
//       }

    }

  }

//   // Assign any remaining nodes to the last domain.
//   for (list<DomainNode<Dimension> >::const_iterator itr = sortedNodeDistribution.begin();
//        itr != sortedNodeDistribution.end();
//        ++itr) {
//     newNodeDistribution.push_back(*itr);
//     newNodeDistribution.back().domainID = numProcs - 1;
//   }

  // Check that all nodes really really have been assigned to a domain.
  BEGIN_CONTRACT_SCOPE;
  CHECK(newNodeDistribution.size() == nodeDistribution.size());
  for (vector<DomainNode<Dimension> >::iterator itr = newNodeDistribution.begin();
       itr != newNodeDistribution.end();
       ++itr) CHECK(itr->domainID >= 0 && itr->domainID < numProcs);
  END_CONTRACT_SCOPE;

  // The nodeDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  CHECK(validDomainDecomposition(newNodeDistribution, dataBase));
  enforceDomainDecomposition(newNodeDistribution, dataBase);

  // Output the final load distribution statistics.
  const string finalLoadStats = this->gatherDomainDistributionStatistics(work);
  if (procID == 0) cerr << "SortAndDivideRedistributeNodes::redistributeNodes final load balancing:" << endl
                        << finalLoadStats << endl << endl;
  MPI_Barrier(Communicator::communicator());

}

//------------------------------------------------------------------------------
// Compute the appropriate number of domains in each dimension for a given
// shape tensor EigenStruct.
//------------------------------------------------------------------------------
vector< vector<int> >
SortAndDivideRedistributeNodes3d::
domainsPerChunk(const Dim<3>::SymTensor::EigenStructType& shapeTensor) const {

  REQUIRE(fuzzyEqual(shapeTensor.eigenValues.magnitude(), 1.0));

  // The total number of domains we are going to assign.
  const int numProcs = this->numDomains();

  // Determine how many x-chunks we're going to use, and the corresponding
  // base number of domains we want in the yz slabs.
  const int xChunks = max(1, int(shapeTensor.eigenValues.x()/shapeTensor.eigenValues.sumElements() * numProcs));
  const int yzChunks = numProcs / xChunks;
  CHECK(xChunks*yzChunks > 0 && xChunks*yzChunks <= numProcs);

  // How many total remainder domains do we have to account for?
  const int totalRemainProcs = numProcs - xChunks*yzChunks;
  CHECK(totalRemainProcs >= 0);

  // Figure out how many remainder domains we have to account for in each
  // yz slab.
  const int baseRemainPerSlab = totalRemainProcs / xChunks;
  vector<int> remainProcs((std::size_t) xChunks, baseRemainPerSlab);
  remainProcs[0] += totalRemainProcs - baseRemainPerSlab * xChunks;
  BEGIN_CONTRACT_SCOPE;
  {
    int checkCount = 0;
    for (vector<int>::const_iterator itr = remainProcs.begin();
         itr != remainProcs.end();
         ++itr) checkCount += *itr;
    CHECK(checkCount == totalRemainProcs);
  }
  END_CONTRACT_SCOPE;

  // Prepare the result.
  vector< vector<int> > result(xChunks);

  // Iterate over each yz slab.
  for (int i = 0; i != xChunks; ++i) {
    CHECK(i < result.size());
    CHECK(i < remainProcs.size());

    // The total number of domains for this slab.
    const int numDomainsInSlab = yzChunks + remainProcs[i];

    // Determine how many y-chunks we're going to use, and the corresponding
    // base number of z-chunks.
    const int yChunks = max(1, int(shapeTensor.eigenValues.y()/(shapeTensor.eigenValues.y() + shapeTensor.eigenValues.z()) * numDomainsInSlab));
    const int zChunks = numDomainsInSlab / yChunks;
    CHECK(yChunks*zChunks > 0 && yChunks*zChunks <= numDomainsInSlab);

    // How many remainder domains do we have to account for?
    const int remainProcs = numDomainsInSlab - yChunks*zChunks;
    CHECK(remainProcs >= 0 && remainProcs < yChunks);

    // Build the numbers of domains in this slab.
    for (int j = 0; j != yChunks; ++j) {
      if (j < remainProcs) {
        result[i].push_back(zChunks + 1);
      } else {
        result[i].push_back(zChunks);
      }
    }
    CHECK(result[i].size() == yChunks);

    BEGIN_CONTRACT_SCOPE;
    {
      int checkCount = 0;
      for (vector<int>::const_iterator itr = result[i].begin();
           itr != result[i].end();
           ++itr) checkCount += *itr;
      CHECK(checkCount == numDomainsInSlab);
    }
    END_CONTRACT_SCOPE;

  }

  // That's it.  Check our post-conditions and return the answer.
  ENSURE(result.size() == xChunks);
  BEGIN_CONTRACT_SCOPE;
  int checkCount = 0;
  for (int i = 0; i != xChunks; ++i) {
    for (vector<int>::const_iterator itr = result[i].begin();
         itr != result[i].end();
         ++itr) checkCount += *itr;
  }
  ENSURE(checkCount == numProcs);
  END_CONTRACT_SCOPE;

  return result;
}

}
}
