//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#include <algorithm>
#include <ctime>

#include "ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/PairComparisons.hh"
namespace Spheral {
namespace NeighborSpace {

using namespace std;

using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using DataBaseSpace::DataBase;
using FieldSpace::FieldList;
using FieldSpace::Field;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConnectivityMap<Dimension>::
ConnectivityMap():
  mNodeLists(),
  mBuildGhostConnectivity(false),
  mConnectivity(),
  mNodeTraversalIndices(),
  mKeys(FieldSpace::FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConnectivityMap<Dimension>::
~ConnectivityMap() {
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
patchConnectivity(const FieldList<Dimension, int>& flags,
                  const FieldList<Dimension, int>& old2new) {

  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();

  // We have to recompute the keys to sort nodes by excluding the 
  // nodes that are being removed.
  const size_t numNodeLists = mNodeLists.size();
  if (domainDecompIndependent) {
    for (size_t iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
      for (size_t i = 0; i != mNodeLists[iNodeList]->numNodes(); ++i) {
        if (flags(iNodeList, i) == 0) mKeys(iNodeList, i) = KeyTraits::maxKey;
      }
    }
  }

  // Iterate over the Connectivity (NodeList).
  vector<size_t> iNodesToKill, jNodesToKill;
  vector<pair<int, Key> > keys, nkeys;
  for (size_t iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
    iNodesToKill = vector<size_t>();
    keys = vector<pair<int, Key> >();

    // Walk the nodes of the NodeList.
    const size_t ioff = mOffsets[iNodeList];
    const size_t numNodes = ((domainDecompIndependent or mBuildGhostConnectivity) ? 
                             mNodeLists[iNodeList]->numNodes() :
                             mNodeLists[iNodeList]->numInternalNodes());

    // Patch the traversal ordering and connectivity for this NodeList.
    for (size_t i = 0; i != numNodes; ++i) {

      // Should we patch this set of neighbors?
      if (flags(iNodeList, i) == 0) {
        iNodesToKill.push_back(i);
      } else {
        if (domainDecompIndependent) keys.push_back(make_pair(old2new(iNodeList, i), mKeys(iNodeList, i)));
        mNodeTraversalIndices[iNodeList][i] = old2new(iNodeList, i);
        vector< vector<int> >& neighbors = mConnectivity[ioff + i];
        CHECK(neighbors.size() == numNodeLists);
        for (size_t jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
          nkeys = vector<pair<int, Key> >();
          jNodesToKill = vector<size_t>();
          for (size_t k = 0; k != neighbors[jNodeList].size(); ++k) {
            const int j = neighbors[jNodeList][k];
            if (flags(jNodeList, j) == 0) {
              jNodesToKill.push_back(k);
            } else {
              if (domainDecompIndependent) nkeys.push_back(make_pair(old2new(jNodeList, j), mKeys(jNodeList, j)));
              neighbors[jNodeList][k] = old2new(jNodeList, j);
            }
          }
          removeElements(neighbors[jNodeList], jNodesToKill);

          // Recompute the ordering of the neighbors.
          if (domainDecompIndependent) {
            sort(nkeys.begin(), nkeys.end(), ComparePairsBySecondElement<pair<int, Key> >());
            for (size_t k = 0; k != neighbors[jNodeList].size(); ++k) {
              CHECK2(k == 0 or nkeys[k].second > nkeys[k-1].second,
                     "Incorrect neighbor ordering:  "
                     << i << " "
                     << k << " "
                     << nkeys[k-1].second << " "
                     << nkeys[k].second);
              neighbors[jNodeList][k] = nkeys[k].first;
            }
          } else {
            sort(neighbors[jNodeList].begin(), neighbors[jNodeList].end());
          }
        }
      }
    }
    removeElements(mNodeTraversalIndices[iNodeList], iNodesToKill);
    // removeElements(*mConnectivity[iNodeList], iNodesToKill);

    // Recompute the ordering for traversing the nodes.
    {
      const size_t numNodes = mNodeTraversalIndices[iNodeList].size();
      if (domainDecompIndependent) {
        // keys = vector<pair<int, Key> >();
        // for (size_t k = 0; k != numNodes; ++k) {
        //   const int i = mNodeTraversalIndices[iNodeList][k];
        //   keys.push_back(make_pair(i, mKeys(iNodeList, i)));
        // }
        sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
        for (size_t k = 0; k != numNodes; ++k) {
          mNodeTraversalIndices[iNodeList][k] = keys[k].first;
        }
      } else {
        for (int i = 0; i != numNodes; ++i) {
          mNodeTraversalIndices[iNodeList][i] = i;
        }
      }
    }
  }
  
  // You can't check valid yet 'cause the NodeLists have not been resized
  // when we call patch!  The valid method should be checked by whoever called
  // this method after that point.
  //ENSURE(valid());
}

//------------------------------------------------------------------------------
// Compute the common neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityIntersectionForNodes(const int nodeListi, const int i,
                                 const int nodeListj, const int j) const {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  const unsigned numNodeLists = mNodeLists.size();
  REQUIRE(nodeListi < numNodeLists and
          nodeListj < numNodeLists);
  const unsigned firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const unsigned firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  REQUIRE(i < firstGhostNodei or j < firstGhostNodej);

  // Prepare the result.
  vector<vector<int> > result(numNodeLists);

  // If both nodes are internal, we simply intersect their neighbor lists.
  if (i < firstGhostNodei and j < firstGhostNodej) {
    vector<vector<int> > neighborsi = this->connectivityForNode(nodeListi, i);
    vector<vector<int> > neighborsj = this->connectivityForNode(nodeListj, j);
    CHECK(neighborsi.size() == numNodeLists);
    CHECK(neighborsj.size() == numNodeLists);
    for (unsigned k = 0; k != numNodeLists; ++k) {
      sort(neighborsi[k].begin(), neighborsi[k].end());
      sort(neighborsj[k].begin(), neighborsj[k].end());
      set_intersection(neighborsi[k].begin(), neighborsi[k].end(),
                       neighborsj[k].begin(), neighborsj[k].end(),
                       back_inserter(result[k]));
    }
  } else if (i < firstGhostNodei) {
    result = this->connectivityForNode(nodeListi, i);
  } else {
    result = this->connectivityForNode(nodeListj, j);
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the union of neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityUnionForNodes(const int nodeListi, const int i,
                          const int nodeListj, const int j) const {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  const unsigned numNodeLists = mNodeLists.size();
  REQUIRE(nodeListi < numNodeLists and
          nodeListj < numNodeLists);
  const unsigned firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const unsigned firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  REQUIRE(i < firstGhostNodei or j < firstGhostNodej);

  // Do the deed.
  vector<vector<int> > result(numNodeLists);
  vector<vector<int> > neighborsi = this->connectivityForNode(nodeListi, i);
  vector<vector<int> > neighborsj = this->connectivityForNode(nodeListj, j);
  CHECK(neighborsi.size() == numNodeLists);
  CHECK(neighborsj.size() == numNodeLists);
  for (unsigned k = 0; k != numNodeLists; ++k) {
    sort(neighborsi[k].begin(), neighborsi[k].end());
    sort(neighborsj[k].begin(), neighborsj[k].end());
    set_union(neighborsi[k].begin(), neighborsi[k].end(),
              neighborsj[k].begin(), neighborsj[k].end(),
              back_inserter(result[k]));
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return the connectivity in terms of global node IDs.
//------------------------------------------------------------------------------
template<typename Dimension>
map<int, vector<int> > 
ConnectivityMap<Dimension>::
globalConnectivity(vector<Boundary<Dimension>*>& boundaries) const {

  // Get the set of global node IDs.
  FieldList<Dimension, int> globalIDs = NodeSpace::globalNodeIDs<Dimension, typename vector<const NodeList<Dimension>*>::const_iterator>
    (mNodeLists.begin(), mNodeLists.end());

  // Make sure all ghost nodes have the appropriate global IDs.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(globalIDs);
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Now convert our connectivity to global IDs.
  map<int, vector<int> > result;
  const size_t numNodeLists = mNodeLists.size();
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    const NodeList<Dimension>* nodeListPtr = mNodeLists[nodeListi];
    for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
      const int gid = globalIDs(nodeListi, i);
      result[gid] = vector<int>();

      const vector< vector<int> >& fullConnectivity = connectivityForNode(nodeListPtr, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];

        for (typename vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) result[gid].push_back(globalIDs(nodeListj, *jItr));

      }
      ENSURE(result[gid].size() == numNeighborsForNode(nodeListPtr, i));
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the index for the given NodeList in our known set.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
ConnectivityMap<Dimension>::
nodeListIndex(const NodeList<Dimension>* nodeListPtr) const {
  return distance(mNodeLists.begin(), 
                  find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConnectivityMap<Dimension>::
valid() const {

  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();

  // Check the offsets.
  const int numNodeLists = mNodeLists.size();
  if (mOffsets.size() != numNodeLists) {
    cerr << "ConnectivityMap::valid: Failed mOffsets.size() == numNodeLists" << endl;
    return false;
  }
  {
    const int numNodes = ((domainDecompIndependent or mBuildGhostConnectivity) ? 
                          mNodeLists.back()->numNodes() : 
                          mNodeLists.back()->numInternalNodes());
    if (mConnectivity.size() != mOffsets.back() + numNodes) {
      cerr << "ConnectivityMap::valid: Failed offset bounding." << endl;
    }
  }

  // Make sure that the NodeLists are listed in the correct sequence, and are
  // FluidNodeLists.
  {
    const NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
    const vector<string> names = registrar.registeredNames();
    int lastPosition = -1;
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      const int newPosition = distance(names.begin(),
                                       find(names.begin(), names.end(), (*itr)->name()));
      if (newPosition <= lastPosition) {
        cerr << "ConnectivityMap::valid: Failed ordering of NodeLists" << endl;
        return false;
      }
      lastPosition = newPosition;
    }
  }

  // Iterate over each NodeList entered.
  int nodeListIDi = 0;
  for (unsigned nodeListIDi = 0; nodeListIDi != numNodeLists; ++nodeListIDi) {

    // Are all internal nodes represented?
    const NodeList<Dimension>* nodeListPtri = mNodeLists[nodeListIDi];
    const int numNodes = ((domainDecompIndependent or mBuildGhostConnectivity) ? 
                          nodeListPtri->numNodes() : 
                          nodeListPtri->numInternalNodes());
    const int firstGhostNodei = nodeListPtri->firstGhostNode();
    if (((nodeListIDi < numNodeLists - 1) and (mOffsets[nodeListIDi + 1] - mOffsets[nodeListIDi] != numNodes)) or
        ((nodeListIDi == numNodeLists - 1) and (mConnectivity.size() - mOffsets[nodeListIDi] != numNodes))) {
      cerr << "ConnectivityMap::valid: Failed test that all nodes set for NodeList "
           << mNodeLists[nodeListIDi]->name()
           << endl;
      return false;
    }

    // Iterate over the nodes for this NodeList.
    const int ioff = mOffsets[nodeListIDi];
    for (int i = 0; i != numNodes; ++i) {

      // The set of neighbors for this node.  This has to be sized as the number of
      // NodeLists.
      const vector< vector<int> >& allNeighborsForNode = mConnectivity[ioff + i];
      if (allNeighborsForNode.size() != numNodeLists) {
        cerr << "ConnectivityMap::valid: Failed allNeighborsForNode.size() == numNodeLists" << endl;
        return false;
      }

      // Iterate over the sets of NodeList neighbors for this node.
      for (int nodeListIDj = 0; nodeListIDj != numNodeLists; ++nodeListIDj) {
        const NodeList<Dimension>* nodeListPtrj = mNodeLists[nodeListIDj];
        const int firstGhostNodej = nodeListPtrj->firstGhostNode();
        const vector<int>& neighbors = allNeighborsForNode[nodeListIDj];

        // We require that the node IDs be sorted, unique, and of course in a valid range.
        if (neighbors.size() > 0) {
          const int minNeighbor = *min_element(neighbors.begin(), neighbors.end());
          const int maxNeighbor = *max_element(neighbors.begin(), neighbors.end());

          if (minNeighbor < 0 or maxNeighbor >= nodeListPtrj->numNodes()) {
            cerr << "ConnectivityMap::valid: Failed test that neighbors must be valid IDs" << endl;
            return false;
          }

          // When enforcing domain independence the ith node may be a ghost, but all of it's neighbors should
          // be internal.
          if (domainDecompIndependent and (i >= firstGhostNodei) and (maxNeighbor > firstGhostNodej)) {
            cerr << "ConnectivityMap::valid: Failed test that all neighbors of a ghost node should be internal." << endl;
            return false;
          }

          for (int k = 1; k < neighbors.size(); ++k) {
            if (domainDecompIndependent) {
              // In the case of domain decomposition reproducibility, neighbors are sorted
              // by hashed IDs.
              if (mKeys(nodeListIDj, neighbors[k]) < mKeys(nodeListIDj, neighbors[k - 1])) {
                cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted for node "
                     << i << endl;
                for (vector<int>::const_iterator itr = neighbors.begin();
                     itr != neighbors.end();
                     ++itr) cerr << *itr << " " 
                                 << mKeys(nodeListIDj, *itr) << " ";
                cerr << endl;
                return false;
              }

            } else {
              // Otherwise they should be sorted by local ID.
              if (neighbors[k] <= neighbors[k - 1]) {
                cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted" << endl;
                for (vector<int>::const_iterator itr = neighbors.begin();
                     itr != neighbors.end();
                     ++itr) cerr << " " << *itr;
                cerr << endl;
                return false;
              }
            }
          }
        }

        // Check that the connectivity is symmetric.
        for (vector<int>::const_iterator jItr = neighbors.begin();
             jItr != neighbors.end();
             ++jItr) {
          if (domainDecompIndependent or mBuildGhostConnectivity or (*jItr < nodeListPtrj->numInternalNodes())) {
            const vector< vector<int> >& otherNeighbors = connectivityForNode(nodeListPtrj, *jItr);
            if (find(otherNeighbors[nodeListIDi].begin(),
                     otherNeighbors[nodeListIDi].end(),
                     i) == otherNeighbors[nodeListIDi].end()) {
              cerr << "ConnectivityMap::valid: Failed test that neighbors must be symmetric." << endl;
              return false;
            }
          }
        }

      }
    }
  }

  // Check that if we are using domain decompostion independence then the keys 
  // have been calculated.
  if (domainDecompIndependent) {
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      if (not mKeys.haveNodeList(**itr)) {
        cerr << "ConnectivityMap::valid: missing information from Keys." << endl;
        return false;
      }
    }
  }

  // Make sure all nodes are listed in the node index traversal stuff.
  if (mNodeTraversalIndices.size() != mNodeLists.size()) {
    cerr << "ConnectivityMap::valid: mNodeTraversalIndices wrong size!" << endl;
    return false;
  }
  for (int nodeList = 0; nodeList != numNodeLists; ++nodeList) {
    const int numExpected = domainDecompIndependent ? mNodeLists[nodeList]->numNodes() : mNodeLists[nodeList]->numInternalNodes();
    bool ok = mNodeTraversalIndices[nodeList].size() == numExpected;
    for (int i = 0; i != numExpected; ++i) {
      ok = ok and (count(mNodeTraversalIndices[nodeList].begin(),
                         mNodeTraversalIndices[nodeList].end(),
                         i) == 1);
    }
    if (not ok) {
      cerr << "ConnectivityMap::valid: mNodeTraversalIndices elements messed up!" << endl;
      return false;
    }
  }

  // Check that the node traversal is ordered correctly.
  for (int nodeList = 0; nodeList != numNodeLists; ++nodeList) {
    if ((domainDecompIndependent and mNodeLists[nodeList]->numNodes() > 0) or
        (not domainDecompIndependent and mNodeLists[nodeList]->numInternalNodes() > 0)) {
      const int firstGhostNode = mNodeLists[nodeList]->firstGhostNode();
      for (const_iterator itr = begin(nodeList); itr < end(nodeList) - 1; ++itr) {
        if (not calculatePairInteraction(nodeList, *itr,
                                         nodeList, *(itr + 1), 
                                         firstGhostNode)) {
          cerr << "ConnectivityMap::valid: mNodeTraversalIndices ordered incorrectly." << endl;
          cerr << *itr << " "
               << *(itr + 1) << " "
               << mKeys(nodeList, *itr) << " "
               << mKeys(nodeList, *(itr + 1)) << " "
               << mNodeLists[nodeList]->positions()(*itr) << " "
               << mNodeLists[nodeList]->positions()(*(itr + 1)) << " "
               << endl;
          for (int i = 0; i != 100; ++i) cerr << mKeys(nodeList, i) << " " << mNodeLists[nodeList]->positions()(i) << " ";
          cerr << endl;
          return false;
        }
      }
    }
  }

  // Everything must be OK.
  return true;
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
computeConnectivity() {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Timing::Time Time;

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      REQUIRE((**itr).neighbor().valid());
    }
    REQUIRE(mOffsets.size() == mNodeLists.size());
  }
  END_CONTRACT_SCOPE

  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  // std::clock_t tpre = std::clock();

  // Build ourselves a temporary DataBase with the set of NodeLists.
  // Simultaneously find the maximum kernel extent.
  DataBase<Dimension> dataBase;
  double kernelExtent = 0.0;
  for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
       itr != mNodeLists.end();
       ++itr) {
    dataBase.appendNodeList(const_cast<NodeList<Dimension>&>(**itr));
    kernelExtent = max(kernelExtent, (**itr).neighbor().kernelExtent());
  }
  const double kernelExtent2 = kernelExtent*kernelExtent;

  // Erase any prior information.
  const unsigned numNodeLists = dataBase.numNodeLists(),
             connectivitySize = mOffsets.back() + 
                                ((domainDecompIndependent or mBuildGhostConnectivity) ? mNodeLists.back()->numNodes() : mNodeLists.back()->numInternalNodes());
  bool ok = (connectivitySize > 0 and mConnectivity.size() == connectivitySize);
  if (ok) {
    CHECK(mNodeTraversalIndices.size() == numNodeLists);
    for (typename ConnectivityStorageType::iterator itr = mConnectivity.begin();
         itr != mConnectivity.end();
         ++itr) {
      CHECK(itr->size() == numNodeLists);
      for (unsigned k = 0; k != numNodeLists; ++k) {
        (*itr)[k].clear();
      }
    }
  } else {
    mConnectivity = ConnectivityStorageType(connectivitySize, vector<vector<int> >(numNodeLists));
    mNodeTraversalIndices = vector<vector<int> >(numNodeLists);
  }

  // If we're trying to be domain decomposition independent, we need a key to sort
  // by that will give us a unique ordering regardless of position.  The Morton ordered
  // hash fills the bill.
  typedef typename KeyTraits::Key Key;
  if (domainDecompIndependent) mKeys = mortonOrderIndices(dataBase);

  // Fill in the ordering for walking the nodes.
  CHECK(mNodeTraversalIndices.size() == numNodeLists);
  if (domainDecompIndependent) {
    for (int iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
      const NodeList<Dimension>& nodeList = *mNodeLists[iNodeList];
      mNodeTraversalIndices[iNodeList].resize(nodeList.numNodes());
      vector<pair<int, Key> > keys;
      keys.reserve(nodeList.numNodes());
      for (int i = 0; i != nodeList.numNodes(); ++i) keys.push_back(pair<int, Key>(i, mKeys(iNodeList, i)));
      sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
      for (int i = 0; i != nodeList.numNodes(); ++i) mNodeTraversalIndices[iNodeList][i] = keys[i].first;
      CHECK(mNodeTraversalIndices[iNodeList].size() == nodeList.numNodes());
    }
  } else {
    for (int iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
      const NodeList<Dimension>& nodeList = *mNodeLists[iNodeList];
      mNodeTraversalIndices[iNodeList].resize(nodeList.numInternalNodes());
      for (int i = 0; i != nodeList.numInternalNodes(); ++i) mNodeTraversalIndices[iNodeList][i] = i;
    }
  }

  // Create a list of flags to keep track of which nodes have been completed thus far.
  FieldList<Dimension, int> flagNodeDone = dataBase.newGlobalFieldList(int());
  flagNodeDone = 0;

  // Get the position and H fields.
  const FieldList<Dimension, Vector> position = dataBase.globalPosition();
  const FieldList<Dimension, SymTensor> H = dataBase.globalHfield();

  // Predeclare stuff we're going to use in the loop.
  unsigned iiNodeList, ii, iNodeList, jNodeList, firstGhostNode;
  int i, firstGhostNodej, j, k, n;
  typename Neighbor<Dimension>::const_iterator masterItr, neighborItr;
  Time start;
  vector<vector<pair<int, Key> > > keys;
  Vector rij;
  Scalar eta2i, eta2j;
  // tpre = std::clock() - tpre;

  // Iterate over the NodeLists.
  // std::clock_t t0, 
  //   tmaster = std::clock_t(0), 
  //   trefine = std::clock_t(0), 
  //   twalk = std::clock_t(0);
  CHECK(mConnectivity.size() == connectivitySize);
  for (iiNodeList = 0; iiNodeList != numNodeLists; ++iiNodeList) {

    // Iterate over the internal nodes in this NodeList, and look
    // for any that are not done yet.
    for (ii = 0; ii != mNodeLists[iiNodeList]->numInternalNodes(); ++ii) {
      if (flagNodeDone(iiNodeList, ii) == 0) {
      
        // Set the master nodes.
        // t0 = std::clock();
        Neighbor<Dimension>::setMasterNeighborGroup(position(iiNodeList, ii),
                                                    H(iiNodeList, ii),
                                                    mNodeLists.begin(),
                                                    mNodeLists.end(),
                                                    mNodeLists[iiNodeList]->neighbor().kernelExtent());
        // tmaster += std::clock() - t0;

        // Iterate over the full of NodeLists again to work on the master nodes.
        for (iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
          const Neighbor<Dimension>& neighbori = mNodeLists[iNodeList]->neighbor();
          CHECK(neighbori.valid());
          Field<Dimension, Scalar>& worki = mNodeLists[iNodeList]->work();
          firstGhostNode = mNodeLists[iNodeList]->firstGhostNode();

          // Iterate over the master nodes in this NodeList.
          for (masterItr = neighbori.masterBegin();
               masterItr != neighbori.masterEnd();
               ++masterItr) {
            i = *masterItr;
            if (domainDecompIndependent or mBuildGhostConnectivity or i < firstGhostNode) {
              CHECK(mOffsets[iNodeList] + i < mConnectivity.size());
              start = Timing::currentTime();
              CHECK(flagNodeDone(iNodeList, i) == 0);

              // Get the neighbor set we're building for this node.
              vector< vector<int> >& neighbors = mConnectivity[mOffsets[iNodeList] + i];
              CHECK2(neighbors.size() == numNodeLists, neighbors.size() << " " << numNodeLists << " " << i);

              // We keep track of the Morton indices.
              keys = vector<vector<pair<int, Key> > >(numNodeLists);

              // Get the state for this node.
              const Vector& ri = position(iNodeList, i);
              const SymTensor& Hi = H(iNodeList, i);

              // Iterate over the neighbor NodeLists.
              for (jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
                Neighbor<Dimension>& neighborj = mNodeLists[jNodeList]->neighbor();
                firstGhostNodej = mNodeLists[jNodeList]->firstGhostNode();

                // Set the refine neighbors.
                // t0 = std::clock();
                neighborj.setRefineNeighborList(ri, Hi);
                // trefine += std::clock() - t0;

                // Iterate over the neighbors in this NodeList.
                // t0 = std::clock();
                for (neighborItr = neighborj.refineNeighborBegin();
                     neighborItr != neighborj.refineNeighborEnd();
                     ++neighborItr) {
                  j = *neighborItr;

                  // Get the neighbor state.
                  const Vector& rj = position(jNodeList, j);
                  const SymTensor& Hj = H(jNodeList, j);

                  // Compute the normalized distance between this pair.
                  rij = ri - rj;
                  eta2i = (Hi*rij).magnitude2();
                  eta2j = (Hj*rij).magnitude2();

                  // If this pair is significant, add it to the list.
                  if (eta2i <= kernelExtent2 or eta2j <= kernelExtent2) {

                    // We don't include self-interactions.
                    if ((iNodeList != jNodeList) or (i != j)) {
                      neighbors[jNodeList].push_back(j);

                      // Do we need ghost connectivity as well?
                      if (j >= firstGhostNodej and (domainDecompIndependent or mBuildGhostConnectivity)) {
                        vector< vector<int> >& otherNeighbors = mConnectivity[mOffsets[jNodeList] + j];
                        CHECK(otherNeighbors.size() == numNodeLists);
                        otherNeighbors[iNodeList].push_back(i);
                      }

                      if (domainDecompIndependent) {
                        keys[jNodeList].push_back(pair<int, Key>(j, mKeys(jNodeList, j)));
                      }

                    }
                  }
                }
                // twalk += std::clock() - t0;
              }
              CHECK(neighbors.size() == numNodeLists);
              CHECK(keys.size() == numNodeLists);
        
              // We have a few options for how to order the neighbors for this node.
              for (k = 0; k != numNodeLists; ++k) {

                if (domainDecompIndependent) {
                  // Sort in a domain independent manner.
                  CHECK(keys[k].size() == neighbors[k].size());
                  sort(keys[k].begin(), keys[k].end(), ComparePairsBySecondElement<pair<int, Key> >());
                  for (j = 0; j != neighbors[k].size(); ++j) neighbors[k][j] = keys[k][j].first;

                } else {
                  // Sort in an attempt to be cache friendly.
                  sort(neighbors[k].begin(), neighbors[k].end());

                }
              }

              // Flag this master node as done.
              flagNodeDone(iNodeList, i) = 1;
              worki(i) += Timing::difference(start, Timing::currentTime());
            }
          }
        }
      }
    }

    // Are we also fleshing out the ghost connectivity?
    if (mBuildGhostConnectivity) {
      iNodeList = iiNodeList;    // Just for consistency.
      firstGhostNode = mNodeLists[iiNodeList]->firstGhostNode();
      n = mNodeLists[iiNodeList]->numNodes();
      for (i = firstGhostNode; i != n; ++i) {
        Neighbor<Dimension>::setMasterNeighborGroup(position(iNodeList, i),
                                                    H(iNodeList, i),
                                                    mNodeLists.begin(),
                                                    mNodeLists.end(),
                                                    mNodeLists[iiNodeList]->neighbor().kernelExtent());
        CHECK(mOffsets[iNodeList] + i < mConnectivity.size());
        CHECK(flagNodeDone(iNodeList, i) == 0);

        // Get the neighbor set we're building for this node.
        vector< vector<int> >& neighbors = mConnectivity[mOffsets[iNodeList] + i];
        CHECK2(neighbors.size() == numNodeLists, neighbors.size() << " " << numNodeLists << " " << i);

        // We keep track of the Morton indices.
        keys = vector<vector<pair<int, Key> > >(numNodeLists);

        // Get the state for this node.
        const Vector& ri = position(iNodeList, i);
        const SymTensor& Hi = H(iNodeList, i);

        // Iterate over the neighbor NodeLists.
        for (jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
          Neighbor<Dimension>& neighborj = mNodeLists[jNodeList]->neighbor();
          firstGhostNodej = mNodeLists[jNodeList]->firstGhostNode();

          // Set the refine neighbors.
          neighborj.setRefineNeighborList(ri, Hi);

          // Iterate over the neighbors in this NodeList.
          for (neighborItr = neighborj.refineNeighborBegin();
               neighborItr != neighborj.refineNeighborEnd();
               ++neighborItr) {
            j = *neighborItr;
            if (j > i) {  // We've already hit all nodes with lower indices than this one.

              // Get the neighbor state.
              const Vector& rj = position(jNodeList, j);
              const SymTensor& Hj = H(jNodeList, j);

              // Compute the normalized distance between this pair.
              rij = ri - rj;
              eta2i = (Hi*rij).magnitude2();
              eta2j = (Hj*rij).magnitude2();

              // If this pair is significant, add it to the list.
              if (eta2i <= kernelExtent2 or eta2j <= kernelExtent2) {

                // We don't include self-interactions.
                if ((iNodeList != jNodeList) or (i != j)) {
                  vector< vector<int> >& otherNeighbors = mConnectivity[mOffsets[jNodeList] + j];
                  CHECK(otherNeighbors.size() == numNodeLists);
                  neighbors[jNodeList].push_back(j);
                  otherNeighbors[iNodeList].push_back(i);
                  if (domainDecompIndependent) {
                    keys[jNodeList].push_back(pair<int, Key>(j, mKeys(jNodeList, j)));
                  }
                }
              }
            }
          }
          CHECK(neighbors.size() == numNodeLists);
          CHECK(keys.size() == numNodeLists);
    
          // We have a few options for how to order the neighbors for this node.
          for (k = 0; k != numNodeLists; ++k) {

            if (domainDecompIndependent) {
              // Sort in a domain independent manner.
              CHECK(keys[k].size() == neighbors[k].size());
              sort(keys[k].begin(), keys[k].end(), ComparePairsBySecondElement<pair<int, Key> >());
              for (j = 0; j != neighbors[k].size(); ++j) neighbors[k][j] = keys[k][j].first;

            } else {
              // Sort in an attempt to be cache friendly.
              sort(neighbors[k].begin(), neighbors[k].end());
            }
          }

          // Flag this master node as done.
          flagNodeDone(iNodeList, i) = 1;
        }
      }
    }
  }

  // In the domain decompostion independent case, we need to sort the neighbors for ghost
  // nodes as well.
  if (domainDecompIndependent) {
    for (int iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
      const NodeList<Dimension>* nodeListPtr = mNodeLists[iNodeList];
      for (int i = nodeListPtr->firstGhostNode();
           i != nodeListPtr->numNodes();
           ++i) {
        vector< vector<int> >& neighbors = mConnectivity[mOffsets[iNodeList] + i];
        CHECK(neighbors.size() == numNodeLists);
        for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
          vector<pair<int, Key> > keys;
          keys.reserve(neighbors[jNodeList].size());
          for (vector<int>::const_iterator itr = neighbors[jNodeList].begin();
               itr != neighbors[jNodeList].end();
               ++itr) keys.push_back(pair<int, Key>(*itr, mKeys(jNodeList, *itr)));
          CHECK(keys.size() == neighbors[jNodeList].size());
          sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
          for (int k = 0; k != keys.size(); ++k) neighbors[jNodeList][k] = keys[k].first;
        }
      }
    }
  }

  // {
  //   tpre = allReduce(unsigned(tpre), MPI_SUM, Communicator::communicator()) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   tmaster = allReduce(unsigned(tmaster), MPI_SUM, Communicator::communicator()) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   trefine = allReduce(unsigned(trefine), MPI_SUM, Communicator::communicator()) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   twalk = allReduce(unsigned(twalk), MPI_SUM, Communicator::communicator()) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   if (Process::getRank() == 0) {
  //     std::cerr << "ConnectivityMap timings (pre, master, refine, walk) = " << tpre << " " << tmaster << " " << trefine << " " << twalk << std::endl;
  //   }
  // }

  // Post conditions.
  BEGIN_CONTRACT_SCOPE
  // Make sure that the correct number of nodes have been completed.
  for (unsigned iNodeList = 0; iNodeList != numNodeLists; ++iNodeList) {
    const unsigned n = (mBuildGhostConnectivity ? 
                        mNodeLists[iNodeList]->numNodes() :
                        mNodeLists[iNodeList]->numInternalNodes());
    for (unsigned i = 0; i != n; ++i) {
      ENSURE2(flagNodeDone(iNodeList, i) == 1,
              "Missed connnectivity for (" << iNodeList << " " << i << ")");
    }
  }
  // Make sure we're ready to be used.
  ENSURE(valid());
  END_CONTRACT_SCOPE
}

}
}

