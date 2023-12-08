#include "computeFragmentField.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/Communicator.hh"

#include <vector>
#include <algorithm>
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
// Reduce a container to it's unique elements.
//------------------------------------------------------------------------------
template<typename ContainerType>
inline
void
reduceToUniqueElements(ContainerType& x) {
  sort(x.begin(), x.end());
  typename ContainerType::iterator itr = unique(x.begin(), x.end());
  x.erase(itr, x.end());
}

//------------------------------------------------------------------------------
// Globally reduce a vector<int> to unique elements across all domains, 
// distributing the result.
//------------------------------------------------------------------------------
inline
void
globalReduceToUniqueElements(vector<int>& x) {

  // Begin by making the local copy unique.
  reduceToUniqueElements(x);

#ifdef USE_MPI
  // If we're parallel, collect the unique set across all processors.
  int procID;
  int numProcs;
  MPI_Comm_rank(Communicator::communicator(), &procID);
  MPI_Comm_size(Communicator::communicator(), &numProcs);
  const vector<int> localX(x);
  x = vector<int>();
  for (int sendID = 0; sendID != numProcs; ++sendID) {
    int n = localX.size();
    MPI_Bcast(&n, 1, MPI_INT, sendID, Communicator::communicator());
    vector<int> otherX;
    if (procID == sendID) {
      otherX = localX;
    } else {
      otherX.resize(n);
    }
    CHECK((int)otherX.size() == n);
    MPI_Bcast(&(*otherX.begin()), n, MPI_INT, sendID, Communicator::communicator());
    x.reserve(x.size() + n);
    copy(otherX.begin(), otherX.end(), back_inserter(x));
  }
  reduceToUniqueElements(x);
  BEGIN_CONTRACT_SCOPE
  {
    int tmp = x.size();
    int sum;
    MPI_Allreduce(&tmp, &sum, 1, MPI_INT, MPI_SUM, Communicator::communicator());
    ENSURE(sum == (int)x.size()*numProcs);
  }
  END_CONTRACT_SCOPE
#endif
}

//------------------------------------------------------------------------------
// computeFragmentField
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, int>
computeFragmentField(const NodeList<Dimension>& nodes,
                     const double linkRadius,
                     const Field<Dimension, typename Dimension::Scalar>& density,
                     const Field<Dimension, typename Dimension::SymTensor>& damage,
                     const double densityThreshold,
                     const double damageThreshold,
                     const bool assignDustToFragments) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  REQUIRE(nodes.numGhostNodes() == 0);
  REQUIRE(density.nodeListPtr() == &nodes);
  REQUIRE(damage.nodeListPtr() == &nodes);

#ifdef USE_MPI
  // Get the rank and total number of processors.
  int procID = 0;
  int numProcs = 1;
  MPI_Comm_rank(Communicator::communicator(), &procID);
  MPI_Comm_size(Communicator::communicator(), &numProcs);

  // Figure out how many elements are in a symmetric tensor.
  SymTensor Hthpt;
  const int Hsize = std::distance(Hthpt.begin(), Hthpt.end());
#endif

  // Grab the state.
  const Field<Dimension, Scalar>& m = nodes.mass();
  const Field<Dimension, Vector>& r = nodes.positions();
  const Field<Dimension, SymTensor>& H = nodes.Hfield();
  Neighbor<Dimension>& neighbor = nodes.neighbor();

  // Get the total number of nodes, and the global IDs on this domain.
  int numGlobalNodesRemaining = numGlobalNodes(nodes);
  vector<int> gIDs;
  {
    Field<Dimension, int> globalNodeField = globalNodeIDs(nodes);
    copy(globalNodeField.begin(), 
         globalNodeField.begin() + nodes.numInternalNodes(),
         back_inserter(gIDs));
  }
  vector<int> globalNodesRemaining(gIDs);
  const vector<int>::iterator maxGlobalItr = max_element(gIDs.begin(), 
                                                         gIDs.end());
  int maxGlobalID = 0;
  if (maxGlobalItr != gIDs.end()) maxGlobalID = *maxGlobalItr;
#ifdef USE_MPI
  {
    int tmp = maxGlobalID;
    MPI_Allreduce(&tmp, &maxGlobalID, 1, MPI_INT, MPI_MAX, Communicator::communicator());
  }
#endif
  maxGlobalID += 1;
  CHECK(maxGlobalID >= numGlobalNodesRemaining);

  // Prepare the result.
  Field<Dimension, int> result(SolidFieldNames::fragmentIDs, nodes, maxGlobalID);
  int numFragments = 0;

  // Flag any nodes above the damage threshold or below the density threshold as dust.
  // Simultaneously remove them from the set of globalNodesRemaining.
  int numDustNodes = 0;
  for (auto i = 0u; i != nodes.numInternalNodes(); ++i) {
    if (damage(i).Trace() > damageThreshold || density(i) < densityThreshold) {
      result(i) = maxGlobalID + 1;
      ++numDustNodes;
      vector<int>::iterator removeItr = find(globalNodesRemaining.begin(),
                                             globalNodesRemaining.end(),
                                             gIDs[i]);
      CHECK(removeItr != globalNodesRemaining.end());
      globalNodesRemaining.erase(removeItr);
    }
  }

  // Reduce the count of remaining nodes by the number of dust nodes.
#ifdef USE_MPI
  {
    int tmp = numDustNodes;
    MPI_Allreduce(&tmp, &numDustNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  }
#endif
  CHECK(numDustNodes >= 0 && numDustNodes <= numGlobalNodesRemaining);
  numGlobalNodesRemaining -= numDustNodes;
  CHECK(numGlobalNodesRemaining >= 0);

  // Reset the neighbor extent for our use.
  const double oldLinkRadius = neighbor.kernelExtent();
  neighbor.kernelExtent(linkRadius);
  neighbor.updateNodes();

  // Iterate until all nodes have been assigned.
  while (numGlobalNodesRemaining > 0) {

    // Find the minimum unassigned node ID.
    const vector<int>::iterator globalMinItr = min_element(globalNodesRemaining.begin(),
                                                           globalNodesRemaining.end());
    int globalMinID = maxGlobalID;
    if (globalMinItr != globalNodesRemaining.end()) globalMinID = *globalMinItr;
#ifdef USE_MPI
    {
      int tmp = globalMinID;
      MPI_Allreduce(&tmp, &globalMinID, 1, MPI_INT, MPI_MIN, Communicator::communicator());
    }
#endif
    CHECK(globalMinID < maxGlobalID);

    // Is this node on this domain?
    int ilocal = globalMinID;
    bool localNode = true;
#ifdef USE_MPI
    int nodeDomain = procID;
    const vector<int>::iterator ilocalItr = find(gIDs.begin(),
                                                 gIDs.end(),
                                                 globalMinID);
    localNode = (ilocalItr != gIDs.end());
    BEGIN_CONTRACT_SCOPE
    {
      int tmp = localNode ? 1 : 0;
      int sum;
      MPI_Allreduce(&tmp, &sum, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CHECK(sum == 1);
    }
    END_CONTRACT_SCOPE
    int tmp = numProcs;
    if (localNode) {
      CHECK(ilocalItr != gIDs.end());
      ilocal = distance(gIDs.begin(), ilocalItr);
      tmp = procID;
      CHECK(result(ilocal) == maxGlobalID);
    }
    MPI_Allreduce(&tmp, &nodeDomain, 1, MPI_INT, MPI_MIN, Communicator::communicator());
    CHECK(nodeDomain >= 0 && nodeDomain < numProcs);
    BEGIN_CONTRACT_SCOPE
    {
      int tmp;
      MPI_Allreduce(&nodeDomain, &tmp, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CHECK(tmp == numProcs*nodeDomain);
    }
    END_CONTRACT_SCOPE
#endif

    // Get the position and H for this node.
    Vector ri;
    SymTensor Hi;
    if (localNode) {
      ri = r(ilocal);
      Hi = H(ilocal);
    }
#ifdef USE_MPI
    MPI_Bcast(&(*ri.begin()), Dimension::nDim, MPI_DOUBLE, nodeDomain, Communicator::communicator());
    MPI_Bcast(&(*Hi.begin()), Hsize, MPI_DOUBLE, nodeDomain, Communicator::communicator());
#endif

    // Find the neighbors for this node within the desired radius.
    vector<int> masterList, coarseNeighbors, refineNeighbors;
    neighbor.setMasterList(ri, Hi, masterList, coarseNeighbors);
    neighbor.setRefineNeighborList(ri, Hi, coarseNeighbors, refineNeighbors);
    vector<int> significantNeighbors;
    vector<int> fragIDs;
    significantNeighbors.reserve(refineNeighbors.size());
    fragIDs.reserve(refineNeighbors.size());
    for (auto itr = refineNeighbors.begin(); itr < refineNeighbors.end(); ++itr) {
      const Vector& rj = r(*itr);
      const SymTensor& Hj = H(*itr);
      const Vector rij = ri - rj;
      const double etai = (Hi*rij).magnitude();
      const double etaj = (Hj*rij).magnitude();
      if (result(*itr) <= maxGlobalID && 
          (etai <= linkRadius || etaj <= linkRadius)) {
        significantNeighbors.push_back(*itr);
        fragIDs.push_back(result(*itr));
      }
    }

    // Distribute the set of fragment ID's.
    globalReduceToUniqueElements(fragIDs);
    CHECK(fragIDs.size() >= 1);

    // Find the minimum fragment ID currently assigned to any of these nodes.
    // If there are no fragments assigned yet, then we'll make this a new fragment ID.
    int fragID = *min_element(fragIDs.begin(), fragIDs.end());
    if (fragID == maxGlobalID) {
      fragID = numFragments;
      numFragments += 1;
    }
    CHECK(fragID >= 0 && fragID < numFragments);
#ifdef USE_MPI
    BEGIN_CONTRACT_SCOPE
    {
      int tmp;
      MPI_Allreduce(&fragID, &tmp, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CHECK(tmp == numProcs*fragID);
      MPI_Allreduce(&numFragments, &tmp, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CHECK(tmp == numProcs*numFragments);
    }
    END_CONTRACT_SCOPE
#endif

    // Remove the known maxGlobalID from the stack of fragment IDs.
    CHECK(fragIDs.back() == maxGlobalID);
    fragIDs.pop_back();

    // Assign the current set of neighbors this fragment ID.
    for (typename vector<int>::const_iterator itr = significantNeighbors.begin();
         itr != significantNeighbors.end();
         ++itr) result(*itr) = fragID;

    // Now assign any nodes that have one of the old fragment IDs to this fragment.
    for (vector<int>::const_iterator fragItr = fragIDs.begin();
         fragItr != fragIDs.end();
         ++fragItr) {
      for (typename Field<Dimension, int>::iterator resultItr = result.begin();
           resultItr != result.end();
           ++resultItr) {
        if (*resultItr == *fragItr) *resultItr = fragID;
      }
    }

    // Make sure we've assigned the node we started with!
    CHECK((localNode && result[ilocal] == fragID) || !localNode);

    // Remove the nodes we've newly assigned from the pool of unassigned global IDs.
    for (typename vector<int>::iterator itr = significantNeighbors.begin();
         itr != significantNeighbors.end();
         ++itr) {
      vector<int>::iterator removeItr = find(globalNodesRemaining.begin(),
                                             globalNodesRemaining.end(),
                                             gIDs[*itr]);
      if (removeItr != globalNodesRemaining.end())
        globalNodesRemaining.erase(removeItr);
    }
    numGlobalNodesRemaining = globalNodesRemaining.size();
#ifdef USE_MPI
    {
      int tmp = numGlobalNodesRemaining;
      MPI_Allreduce(&tmp, &numGlobalNodesRemaining, 1, MPI_INT, MPI_SUM, Communicator::communicator());
    }
#endif

    BEGIN_CONTRACT_SCOPE
    {
      for (typename vector<int>::iterator itr = significantNeighbors.begin();
           itr != significantNeighbors.end();
           ++itr) {
        CHECK(find(globalNodesRemaining.begin(),
                   globalNodesRemaining.end(),
                   gIDs[*itr]) == globalNodesRemaining.end());
      }
    }
    END_CONTRACT_SCOPE

  }

  // Make sure all nodes have been assigned to a valid fragment.
  CHECK(nodes.numInternalNodes() == 0 ||
        *min_element(result.begin(), result.end()) >= 0);

//   // Assign the dust nodes a fragment index of -1.
//   for (int i = 0; i != nodes.numInternalNodes(); ++i) {
//     if (result(i) == maxGlobalID + 1) result(i) = -1;
//   }

  // At this point all nodes have been assigned to unique fragment IDs, but those IDs are not
  // necessarily contiguous.  Go ahead and make them a contiguous set.
  vector<int> fragIDs(result.begin(), result.end());
  globalReduceToUniqueElements(fragIDs);
  numFragments = fragIDs.size();
  for (int i = 0; i != numFragments; ++i) {
    for (auto j = 0u; j != nodes.numInternalNodes(); ++j) {
      if (result(j) == fragIDs[i]) result(j) = i;
    }
  }

  // If requested, assign the dust to the nearest fragment (as defined by center
  // of mass distances).
  if (assignDustToFragments) {
    const int dustID = numFragments - 1;

    // Find the center of mass for each fragment.
    vector<Vector> rfrag(numFragments - 1);
    vector<double> mfrag(numFragments - 1);
    for (auto i = 0u; i != nodes.numInternalNodes(); ++i) {
      if (result[i] != dustID) {
        CHECK(distinctlyGreaterThan(m(i), 0.0));
        mfrag[result[i]] += m(i);
        rfrag[result[i]] += m(i)*r(i);
      }
    }
#ifdef USE_MPI
    for (int i = 0; i != numFragments - 1; ++i) {
      double mtmp = mfrag[i];
      Vector rtmp = rfrag[i];
      MPI_Allreduce(&mtmp, &(mfrag[i]), 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
      MPI_Allreduce(&(*rtmp.begin()), &(*rfrag[i].begin()), Dimension::nDim, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
    }
#endif
    for (int i = 0; i != numFragments - 1; ++i) {
      CHECK(distinctlyGreaterThan(mfrag[i], 0.0));
      rfrag[i] /= mfrag[i];
    }

    // Now go over all the dust nodes, find the fragment they're closeet to,
    // and assign them to that fragment.
    for (auto i = 0u; i != nodes.numInternalNodes(); ++i) {
      if (result[i] == dustID) {
        const Vector& ri = r(i);
        double rmin = DBL_MAX;
        if (numFragments > 1) {
          int fragmin = -1;
          for (int j = 0; j != numFragments - 1; ++j) {
            const double dr2 = (ri - rfrag[j]).magnitude2();
            if (dr2 < rmin) {
              rmin = dr2;
              fragmin = j;
            }
          }
          CHECK(fragmin >= 0 && fragmin < numFragments - 1);
          result[i] = fragmin;
        } else {
          result[i] = 0;
        }
      }
    }

    // Check that all dust nodes have been assigned.
    for (auto i = 0u; i != nodes.numInternalNodes(); ++i) 
      CHECK(numFragments == 1 or result[i] < dustID);
    
  }

  // Set the neighbor state back how we found it.
  neighbor.kernelExtent(oldLinkRadius);
  neighbor.updateNodes();

  return result;
}

}
