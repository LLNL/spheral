//---------------------------------Spheral++----------------------------------//
// SpaceFillingCurveRedistributeNodes
//
// This is an abstract base for the space filling curve family of 
// repartitioners.  The assumption is that the descendent classes will provide
// the computeHashedIndices method to assign unique keys to each point in the
// order that that algorithm wants the points distributed.
//
// Created by JMO, Wed Apr  9 13:13:46 PDT 2008
//----------------------------------------------------------------------------//
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <bitset>

#include "SpaceFillingCurveRedistributeNodes.hh"
#include "DomainNode.hh"
#include "DistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/RedistributionRegistrar.hh"
#include "Utilities/allReduce.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace PartitionSpace {

using namespace std;

using DataBaseSpace::DataBase;
using NodeSpace::NodeList;
using BoundarySpace::DistributedBoundary;
using BoundarySpace::Boundary;
using FieldSpace::FieldList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Compare pairs of <Key, *> by the first element (Key).
//------------------------------------------------------------------------------
template<typename DataType>
struct SortNodesByHashedIndex {
  bool operator()(const pair<unsigned long long, DataType>& lhs,
                  const pair<unsigned long long, DataType>& rhs) {
    return lhs.first < rhs.first;
  }
};

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpaceFillingCurveRedistributeNodes<Dimension>::
SpaceFillingCurveRedistributeNodes(double dummy,
                                   const double minNodesPerDomainFraction,
                                   const double maxNodesPerDomainFraction,
                                   const bool workBalance,
                                   const bool localReorderOnly):
  RedistributeNodes<Dimension>(),
  mMinNodesPerDomainFraction(minNodesPerDomainFraction),
  mMaxNodesPerDomainFraction(maxNodesPerDomainFraction),
  mWorkBalance(workBalance),
  mLocalReorderOnly(localReorderOnly) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SpaceFillingCurveRedistributeNodes<Dimension>::
~SpaceFillingCurveRedistributeNodes() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Call on the SpaceFillingCurve library to describe an
// optimal partioning of the nodes, and then apply that partitioning.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> boundaries) {

  // The usual parallel info.
  const int numProcs = this->numDomains();
  const int procID = this->domainID();

  // The min and max range of nodes we want per domain.
  const int totalNumNodes = dataBase.globalNumInternalNodes();
  const int minNodes = int(mMinNodesPerDomainFraction * double(totalNumNodes)/numProcs);
  const int maxNodes = int(mMaxNodesPerDomainFraction * double(totalNumNodes)/numProcs);

  // Below a certain level of work per process we just punt.
  if (double(totalNumNodes)/double(Process::getTotalNumberOfProcesses()) >= 1.0) {

    // Get the global IDs.
    const FieldList<Dimension, int> globalIDs = NodeSpace::globalNodeIDs(dataBase);

    // Compute the work per node.
    FieldList<Dimension, Scalar> workField(FieldSpace::Copy);
    if (this->workBalance()) {

      // Enforce boundary conditions for the work computation.
      for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
           nodeListItr != dataBase.nodeListEnd();
           ++nodeListItr) {
        (*nodeListItr)->numGhostNodes(0);
        (*nodeListItr)->neighbor().updateNodes();
      }
      for (typename vector<Boundary<Dimension>*>::iterator boundaryItr = boundaries.begin(); 
           boundaryItr != boundaries.end();
           ++boundaryItr) {
        (*boundaryItr)->setAllGhostNodes(dataBase);
        (*boundaryItr)->finalizeGhostBoundary();
        for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
             nodeListItr != dataBase.fluidNodeListEnd(); 
             ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
      }

      // Get the local description of the domain distribution, with the work per node filled in.
      workField = this->workPerNode(dataBase, 1.0);

    } else {

      // We assign a constant work of 1 to each node, effectively balancing by
      // node count.
      workField = dataBase.newGlobalFieldList(1.0, "work");

    }

    // Print the beginning statistics.
    std::string stats0 = this->gatherDomainDistributionStatistics(workField);
    if (Process::getRank() == 0) cout << "SpaceFillingCurveRedistributeNodes: INITIAL node distribution statistics:" << endl
                                      << stats0 << endl;

    // Now we can get the node distribution description.
    vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs, workField);

    // Clear out any ghost nodes.
    // We won't bother to update the neighbor info at this point -- we don't need 
    // it for this algorithm, so we just update it when we're done.
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
    }

    // Compute the target work per domain.
    const Scalar targetWork = workField.sumElements()/numProcs;
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: Target work per process " << targetWork << endl;

    // Compute the Key indices for each point on this processor.
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: Hashing indices" << endl;
    FieldList<Dimension, Key> indices = computeHashedIndices(dataBase);

    // Find the range of hashed indices.
    const Key indexMin = indices.min();
    const Key indexMax = indices.max();
    CHECK(indexMax < indexMax + indexMax);
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: Index min/max : " << indexMin << " " << indexMax << endl;

    // Build the array of (hashed index, DomainNode) pairs.
    // Note this comes back locally sorted.
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: sorting indices" << endl;
    vector<pair<Key, DomainNode<Dimension> > > sortedIndices = buildIndex2IDPairs(indices,
                                                                                  nodeDistribution);
    const int numLocalNodes = nodeDistribution.size();

    // Build our set of unique indices and their count.
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: Counting uniques and such" << endl;
    vector<Key> uniqueIndices;
    vector<int> count;
    vector<Scalar> work;
    int maxCount = 1;
    if (sortedIndices.size() > 0) {
      uniqueIndices.reserve(sortedIndices.size());
      count.reserve(sortedIndices.size());
      work.reserve(sortedIndices.size());
      uniqueIndices.push_back(sortedIndices.front().first);
      count.push_back(1);
      work.push_back(sortedIndices.front().second.work);
      {
        int j = 0;
        for (int i = 1; i < sortedIndices.size(); ++i) {
          CHECK(sortedIndices[i].first >= sortedIndices[i-1].first);
          if (sortedIndices[i].first == uniqueIndices[j]) {
            ++count[j];
            work[j] += sortedIndices[i].second.work;
          } else {
            uniqueIndices.push_back(sortedIndices[i].first);
            count.push_back(1);
            work.push_back(sortedIndices[i].second.work);
            ++j;
          }
          maxCount = max(maxCount, count[j]);
        }
      }
      CHECK(count.size() == uniqueIndices.size());
      CHECK(work.size() == uniqueIndices.size());
    }
    maxCount = allReduce(maxCount, MPI_MAX, Communicator::communicator());
    if (procID == 0) cerr << "SpaceFillingCurveRedistributeNodes: max redundancy is " << maxCount << endl;

    //   // DEBUG
    //   {
    //     for (int ii = 0; ii != numProcs; ++ii) {
    //       if (procID == ii) {
    //         for (size_t i = 0; i != uniqueIndices.size(); ++i) {
    //           cerr << uniqueIndices[i] << endl;
    //         }
    //       }
    //       MPI_Barrier(Communicator::communicator());
    //     }
    //   }
    //   // DEBUG

    // Figure out the range of hashed indices we want for each process.
    // Note this will not be optimal when there are degnerate indices!
    Key lowerBound = indexMin;
    vector<pair<Key, Key> > indexRanges;
    for (int iProc = 0; iProc != numProcs; ++iProc) {
      Key upperBound;
      int numNodes;
      findUpperKey(uniqueIndices,
                   count,
                   work,
                   lowerBound,
                   indexMax,
                   targetWork,
                   minNodes,
                   maxNodes,
                   upperBound,
                   numNodes);
      // cerr << "Chose indices in range: "
      //      << lowerBound << " " << upperBound << endl;
      // if (procID == 0) cerr << "  range [" 
      //                       << lowerBound << " " 
      //                       << upperBound << "], work in range "
      //                       << workInRange(uniqueIndices, work, lowerBound, upperBound)
      //                       << ", num in range "
      //                       << numIndicesInRange(uniqueIndices, count, lowerBound, upperBound)
      //                       << endl;
      indexRanges.push_back(pair<Key, Key>(lowerBound, upperBound));
      lowerBound = findNextIndex(uniqueIndices, upperBound, indexMax);
    }
    CHECK(lowerBound == indexMax);
    CHECK(indexRanges[0].first == indexMin);
    CHECK(indexRanges.back().second == indexMax);
    CHECK(indexRanges.size() == numProcs);

    // We now know the target index range for each domain.
    // Go through our local DomainNode set and assign them appropriately.
    nodeDistribution = vector<DomainNode<Dimension> >();
    for (typename vector<pair<Key, DomainNode<Dimension> > >::iterator itr = sortedIndices.begin();
         itr != sortedIndices.end();
         ++itr) {
      nodeDistribution.push_back(itr->second);
      nodeDistribution.back().domainID = domainForIndex(itr->first, indexRanges);
      CHECK(nodeDistribution.back().domainID >= 0 and
            nodeDistribution.back().domainID < numProcs);
    }
    CHECK(nodeDistribution.size() == numLocalNodes);

    // Redistribute nodes between domains.
    CHECK(this->validDomainDecomposition(nodeDistribution, dataBase));
    if (not this->localReorderOnly()) this->enforceDomainDecomposition(nodeDistribution, dataBase);

    // At this point we have the nodes on each domain, but they are not sorted locally.
    // That is the next step.  Rebuild the local sorting.
    nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs);
    sortedIndices = buildIndex2IDPairs(indices, nodeDistribution);

    // Extract the desired node orderings for each NodeList.
    const size_t numNodeLists = dataBase.numNodeLists();
    vector<vector<int> > orderings;
    vector<int> newNodeIDs(numNodeLists, 0);
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) orderings.push_back(vector<int>((*nodeListItr)->numInternalNodes()));
    for (typename vector<pair<Key, DomainNode<Dimension> > >::const_iterator itr = sortedIndices.begin();
         itr != sortedIndices.end();
         ++itr) {
      const DomainNode<Dimension>& node = itr->second;
      CHECK(node.nodeListID >= 0 and node.nodeListID < numNodeLists);
      CHECK(node.localNodeID >= 0 and node.localNodeID < orderings[node.nodeListID].size());
      CHECK(newNodeIDs[node.nodeListID] < orderings[node.nodeListID].size());
      orderings[node.nodeListID][newNodeIDs[node.nodeListID]] = node.localNodeID;
      ++newNodeIDs[node.nodeListID];
    }

    // Now we have the local orderings, so enforce them!
    {
      int nodeListID = 0;
      for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
           nodeListItr != dataBase.nodeListEnd();
           ++nodeListItr, ++nodeListID) (*nodeListItr)->reorderNodes(orderings[nodeListID]);
    }

    // Reinitialize neighbor info.
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) {
      (*nodeListItr)->neighbor().updateNodes();
    }

    // Notify everyone that the nodes have just been shuffled around.
    RedistributionRegistrar::instance().broadcastRedistributionNotifications();

    // Print the final statistics.
    std::string stats1 = this->gatherDomainDistributionStatistics(workField);
    if (Process::getRank() == 0) cout << "SpaceFillingCurveRedistributeNodes: FINAL node distribution statistics:" << endl
                                      << stats1 << endl;

    // Post-conditions.
    // Make sure the nodes are now sorted by the index keys.
    BEGIN_CONTRACT_SCOPE
      {
        for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
             nodeListItr != dataBase.nodeListEnd();
             ++nodeListItr) {
          const Field<Dimension, Key> keyField = **indices.fieldForNodeList(**nodeListItr);
          for (int i = 1; i < (*nodeListItr)->numInternalNodes(); ++i) {
            ENSURE(keyField(i) >= keyField(i - 1));
          }
        }
      }
    END_CONTRACT_SCOPE
  }
}


//------------------------------------------------------------------------------
// Based on the bounding box and our hardwired number of cells, compute the
// step size in each dimension.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
SpaceFillingCurveRedistributeNodes<Dimension>::
computeStepSize(const pair<typename Dimension::Vector, typename Dimension::Vector>& box) const {
  // We choose the number of cells per dimension to be the maximum power of 2 such that
  // mNumGridCellsPerSide**3 < 2**64, so
  // that the indices will fit into a 64 bit Key.  This gives us 2**21.
  const Vector result = (box.second - box.first)/1048576;
//   const Vector result = (box.second - box.first)/2097152;
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Stitch together the set of indices and DomainNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<pair<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key, DomainNode<Dimension> > >
SpaceFillingCurveRedistributeNodes<Dimension>::
buildIndex2IDPairs(const FieldList<Dimension, typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                   const vector<DomainNode<Dimension> >& domainNodes) const {

  vector<pair<Key, DomainNode<Dimension> > > result;
  result.reserve(domainNodes.size());
  size_t i = 0;
  for (InternalNodeIterator<Dimension> nodeItr = indices.internalNodeBegin();
       nodeItr != indices.internalNodeEnd();
       ++nodeItr, ++i) {
    CHECK(i < domainNodes.size());
    result.push_back(pair<Key, DomainNode<Dimension> >(indices(nodeItr),
                                                                 domainNodes[i]));
  }
  ENSURE(i == domainNodes.size());
  ENSURE(result.size() == domainNodes.size());

  // Sort locally by the hashed index.
  sort(result.begin(), result.end(), SortNodesByHashedIndex<DomainNode<Dimension> >());
  return result;
}

//------------------------------------------------------------------------------
// Find the (global) hashed index the given amount of work above the specified
// lower bound.
// We assume here that the vector<indices> is unique, and the associated 
// work array contains the work per node for each index.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SpaceFillingCurveRedistributeNodes<Dimension>::Key
SpaceFillingCurveRedistributeNodes<Dimension>::
findUpperKey(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
             const vector<int>& count,
             const vector<typename Dimension::Scalar>& work,
             const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
             const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key maxUpperBound,
             const typename Dimension::Scalar workTarget,
             const int minNodes,
             const int maxNodes,
             Key& upperBound1,
             int& numNodes) const {
  REQUIRE(count.size() == indices.size());
  REQUIRE(work.size() == indices.size());
  REQUIRE(lowerBound <= maxUpperBound);
  REQUIRE(workTarget > 0.0);
  
  // Now bisect in on the upper bound.
  Key upperBound0 = lowerBound;
  upperBound1 = maxUpperBound;
  while (upperBound1 - upperBound0 > KeyTraits::one) {
    const Key trialUpperBound = (upperBound0 + upperBound1)/KeyTraits::two;
    CHECK(trialUpperBound >= upperBound0 and
          trialUpperBound <= upperBound1);
    Scalar workTrial;
    workAndNodesInRange(indices, count, work, lowerBound, trialUpperBound, numNodes, workTrial);
    // if (Process::getRank() == 0) cerr << "  --> [" 
    //                                   << upperBound0 << " " << trialUpperBound 
    //                                   << "]   "
    //                                   << workTrial << " " << workTarget 
    //                                   << endl;
    if (workTrial < workTarget) {
      if (numNodes <= maxNodes) {
        upperBound0 = trialUpperBound;
      } else {
        upperBound1 = trialUpperBound;
      }
    } else {
      if (numNodes >= minNodes) {
        upperBound1 = trialUpperBound;
      } else {
        upperBound0 = trialUpperBound;
      }
    }
  }

  // That's it.
  ENSURE(upperBound1 <= maxUpperBound);
  return upperBound1;
}

//------------------------------------------------------------------------------
// Compute the (global) number of nodes in the given range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
numIndicesInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                   const vector<int>& count,
                   const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
                   const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound) const {
  REQUIRE(lowerBound <= upperBound);

  int result = 0;
  if (indices.size() > 0) {

    // Find the positions of the requested range in our local set of indices.
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, std::min(int(indices.size()) - 1, bisectSearch(indices, upperBound)));

    // Count up how many of the suckers we have locally.
    for (int i = ilower; i != iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) result += count[i];
    }
  }

  // Globally reduce that sucker.
  result = allReduce(result, MPI_SUM, Communicator::communicator());

  return result;
}

//------------------------------------------------------------------------------
// Compute the (global) work for the nodes in the given range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
SpaceFillingCurveRedistributeNodes<Dimension>::
workInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
            const vector<typename Dimension::Scalar>& work,
            const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
            const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound) const {
  REQUIRE(lowerBound <= upperBound);

  Scalar result = 0;

  if (indices.size() > 0) {
    // Find the positions of the requested range in our local set of indices.
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, min(int(int(indices.size()) - 1), bisectSearch(indices, upperBound)));

    // Sum our local work.
    for (int i = ilower; i != iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) result += work[i];
    }
  }

  // Globally reduce that sucker.
  result = allReduce(result, MPI_SUM, Communicator::communicator());

  return result;
}

//------------------------------------------------------------------------------
// Compute the (global) work and number of nodes for the nodes in the given
//  range of indices.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
workAndNodesInRange(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
                    const vector<int>& count,
                    const vector<typename Dimension::Scalar>& work,
                    const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key lowerBound,
                    const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key upperBound,
                    int& countInRange,
                    Scalar& workInRange) const {

  REQUIRE(lowerBound <= upperBound);

  // Find the positions of the requested range in our local set of indices.
  workInRange = 0.0;
  countInRange = 0;
  if (indices.size() > 0) {
    const int ilower = max(0, bisectSearch(indices, lowerBound));
    const int iupper = max(0, min(int(int(indices.size()) - 1), bisectSearch(indices, upperBound)));

    // Sum our local work.
    for (int i = ilower; i < iupper + 1; ++i) {
      if (indices[i] >= lowerBound and 
          indices[i] <= upperBound) {
        workInRange += work[i];
        countInRange += count[i];
      }
    }
  }

  // Globally reduce that sucker.
  workInRange = allReduce(workInRange, MPI_SUM, Communicator::communicator());
  countInRange = allReduce(countInRange, MPI_SUM, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Compute the number of nodes we want on the given domain.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
targetNumNodes(const int numGlobal,
               const int numProcs,
               const int targetProc) const {
  const int nlocalp = numGlobal/numProcs;
  const int nremainder = numGlobal % numProcs;
  const int result = targetProc > nremainder ? nlocalp : (nlocalp + 1);
  ENSURE(result > 0 and result <= numGlobal);
  return result;
}

//------------------------------------------------------------------------------
// Find the next unique index globally following the given one.
// Note that if there is no next index, we return the maximum value.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SpaceFillingCurveRedistributeNodes<Dimension>::Key
SpaceFillingCurveRedistributeNodes<Dimension>::
findNextIndex(const vector<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key>& indices,
              const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key index,
              const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key maxIndex) const {

  // Find the position of the given index in our local array.
  const int inext = bisectSearch(indices, index) + 1;
  CHECK(inext >= 0 and inext <= indices.size());

  // The next value on this processor.
  Key result = maxIndex;
  if (indices.size() > 0) {
    if (inext < indices.size() and indices[inext] > index) result = indices[inext];
    //   cerr << "  Local result: " << result << endl;

    // Get the global answer.
    result = allReduce(result, MPI_MIN, Communicator::communicator());
    //   cerr << "Global result: " << result << endl;
  } else {
    result = 0;
  }

  // That's it.
  ENSURE(result <= maxIndex);
  return result;
}

//------------------------------------------------------------------------------
// Find the domain for the given index given the set of index ranges for 
// each processors.
//------------------------------------------------------------------------------
template<typename Dimension>
int
SpaceFillingCurveRedistributeNodes<Dimension>::
domainForIndex(const typename SpaceFillingCurveRedistributeNodes<Dimension>::Key index,
               const vector<pair<typename SpaceFillingCurveRedistributeNodes<Dimension>::Key, 
                                 typename SpaceFillingCurveRedistributeNodes<Dimension>::Key> >& indexRanges) const {

  // Find that sucker by bisection.
  const int numDomains = indexRanges.size();
  int ilower = 0;
  int iupper = numDomains;
  while (iupper - ilower > 1) {
    const int imid = (iupper + ilower)/2;
    if (indexRanges[imid].first > index) {
      iupper = imid;
    } else {
      ilower = imid;
    }
  }

  // Thats it.
  ENSURE(indexRanges[ilower].first <= index and
         indexRanges[ilower].second >= index);
  return ilower;
}

//------------------------------------------------------------------------------
// Range of allowed nodes per domain as a fraction of the average number per 
// domain.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SpaceFillingCurveRedistributeNodes<Dimension>::
minNodesPerDomainFraction() const {
  return mMinNodesPerDomainFraction;
}

template<typename Dimension>
double
SpaceFillingCurveRedistributeNodes<Dimension>::
maxNodesPerDomainFraction() const {
  return mMaxNodesPerDomainFraction;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
minNodesPerDomainFraction(const double x) {
  mMinNodesPerDomainFraction = x;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
maxNodesPerDomainFraction(const double x) {
  mMaxNodesPerDomainFraction = x;
}

//------------------------------------------------------------------------------
// Flag for whether we should compute the work per node or strictly balance by
// node count.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpaceFillingCurveRedistributeNodes<Dimension>::
workBalance() const {
  return mWorkBalance;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
workBalance(const bool val) {
  mWorkBalance = val;
}

//------------------------------------------------------------------------------
// Flag that will cause us not to repartition between domains, but only sort locally on 
// each domain.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpaceFillingCurveRedistributeNodes<Dimension>::
localReorderOnly() const {
  return mLocalReorderOnly;
}

template<typename Dimension>
void
SpaceFillingCurveRedistributeNodes<Dimension>::
localReorderOnly(const bool val) {
  mLocalReorderOnly = val;
}

}
}

