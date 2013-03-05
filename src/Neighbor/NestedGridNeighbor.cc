//---------------------------------Spheral++----------------------------------//
// NestedGridNeighbor -- A reimplementation of the nested AMR like neighboring 
// method used in the original fortran version of Spheral.  Described in
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by JMO, Wed Dec 22 21:02:47 PST 1999
//----------------------------------------------------------------------------//
#include <math.h>
#include <limits.h>
#include <float.h>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

#include "NestedGridNeighbor.hh"
#include "Neighbor.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "GridCellIndex.hh"
#include "Boundary/mapPositionThroughPlanes.hh"
#include "Utilities/intpow2.hh"
#include "DBC.hh"

namespace Spheral {
namespace NeighborSpace {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with everything necessary to completely specify the neighboring.
//------------------------------------------------------------------------------
template<typename Dimension>
NestedGridNeighbor<Dimension>::
NestedGridNeighbor(NodeList<Dimension>& aNodeList,
                   const NeighborSearchType aSearchType,
                   const int numGridLevels,
                   const double topGridCellSize,
                   const typename Dimension::Vector origin,
                   const double kernelExtent,
		   const int gridCellInfluenceRadius):
  Neighbor<Dimension>(aNodeList, aSearchType, kernelExtent),
  mMaxGridLevels(numGridLevels),
  mFirstParentGridLevel(0),
  mGridCellInfluenceRadius(gridCellInfluenceRadius),
  mGridOrigin(), // (origin),
  mGridLevelOccupied(numGridLevels, 0),
  mGridLevelConst0(log(topGridCellSize * gridCellInfluenceRadius)*ln2inverse),
  mGridCellSizeInv(numGridLevels),
  mGridCellHead(numGridLevels),
  mNodeInCell(numGridLevels),
  mNextNodeInCell(aNodeList.numNodes(), mEndOfLinkList),
  mNodeOnGridLevel(aNodeList.numNodes()),
  mMasterGridLevel(-1),
  mMasterGridCellIndex(),
//   mDaughterCells(numGridLevels),
//   mEmptyNest(),
  mOccupiedGridCells(numGridLevels) {

  // CHECK the input is sane.
  VERIFY(numGridLevels > 0 and numGridLevels < 32);
  VERIFY(topGridCellSize > 0.0);
  VERIFY(kernelExtent > 0.0);
  
  this->kernelExtent(kernelExtent);

  // Build the vector of gridcell sizes.
  for (int i = 0; i < numGridLevels; ++i) mGridCellSizeInv[i] = double(intpow2(i))/topGridCellSize;

  // We should have enough info, so go ahead and build the node info.
  updateNodes();

  // When this constructor is through, the neighbor info should be ready to use.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NestedGridNeighbor<Dimension>::
~NestedGridNeighbor() {
}

//------------------------------------------------------------------------------
// Set the master list of nodes, required by the base classes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setMasterList(const typename Dimension::Vector& position,
              const typename Dimension::Scalar& H) {
  setNestedMasterList(position, H);
}

template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setMasterList(const typename Dimension::Vector& position,
              const typename Dimension::SymTensor& H) {
  setNestedMasterList(position, H);
}

template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setMasterList(const typename Dimension::Vector& position) {
  SymTensor fakeH;
  fakeH.Identity();
  fakeH *= 1.0e30;
  setNestedMasterList(position, fakeH);
}

//------------------------------------------------------------------------------
// Set the refined list of nodes, required by the base classes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setRefineNeighborList(const typename Dimension::Vector& position,
                      const typename Dimension::Scalar& H) {
  setNestedRefineNeighborList(position, H);
}

template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setRefineNeighborList(const typename Dimension::Vector& position,
                      const typename Dimension::SymTensor& H) {
  setNestedRefineNeighborList(position, H);
}

template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setRefineNeighborList(const typename Dimension::Vector& position) {
  SymTensor fakeH;
  fakeH.Identity();
  fakeH *= 1.0e30;
  setNestedRefineNeighborList(position, fakeH);
}

//------------------------------------------------------------------------------
// Set the master list of nodes based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setMasterList(const int nodeID) {
  REQUIRE(valid());
  REQUIRE(nodeID >= 0 and nodeID < this->nodeList().numInternalNodes());
  vector<int>& masterList = this->accessMasterList();
  masterList = vector<int>();
  const Field<Dimension, Vector>& positions = this->nodeList().positions();
  const Field<Dimension, SymTensor>& Hfield = this->nodeList().Hfield();
  BEGIN_CONTRACT_SCOPE;
  {
    const int gridLevelCheck = mNodeOnGridLevel[nodeID];
    const GridCellIndex<Dimension> gridCellCheck = mNodeInCell[gridLevelCheck][nodeID];
    CHECK(gridLevel(Hfield(nodeID)) == gridLevelCheck);
    CHECK(gridCellIndex(positions(nodeID), gridLevel(nodeID)) == gridCellCheck);
  }
  END_CONTRACT_SCOPE;
  setMasterList(positions(nodeID), Hfield(nodeID));
}

//------------------------------------------------------------------------------
// Set the refined list of potential neighbors based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setRefineNeighborList(const int nodeID) {
  REQUIRE(valid());
  REQUIRE(nodeID >= 0 and nodeID < this->nodeList().numInternalNodes());
  REQUIRE(find(this->accessMasterList().begin(), this->accessMasterList().end(), nodeID) !=
        this->accessMasterList().end());
  const Field<Dimension, Vector>& positions = this->nodeList().positions();
  const Field<Dimension, SymTensor>& Hfield = this->nodeList().Hfield();
  setRefineNeighborList(positions(nodeID), Hfield(nodeID));
}

//------------------------------------------------------------------------------
// Set the master list of nodes for which we can quickly compute neighbors.
// For NestedGridNeighbor objects, this is just all the nodes that belong to
// the given node (position and extent).
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename OtherHType>
void
NestedGridNeighbor<Dimension>::
setNestedMasterList(const typename Dimension::Vector& position,
                    const OtherHType& H) {
  const int gl = gridLevel(H);
  const GridCellIndex<Dimension> gc = gridCellIndex(position, gl);
  setNestedMasterList(gc, gl);
}

//------------------------------------------------------------------------------
// Set the master list of nodes for which we can quickly compute neighbors.
// For NestedGridNeighbor objects, this is just all the nodes that belong to
// the given node's master grid cell.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setNestedMasterList(const GridCellIndex<Dimension>& gridCell,
                    const int gridLevel) {

  REQUIRE(valid());

  // Set the gridlevel and gridcell.
  mMasterGridCellIndex = gridCell;
  mMasterGridLevel = gridLevel;

  // Set the master list to all (internal) nodes in this gridcell (on this gridlevel).
  vector<int>& mMasterList = this->accessMasterList();
  mMasterList = internalNodesInCell(mMasterGridCellIndex, mMasterGridLevel);

  // We also need to set the coarse neighbor list -- potential neighbors for all
  // nodes in the master gridcell.  We can outsource this to our private findNestedNeighbors
  // method.
  vector<int>& coarseNeighborList = this->accessCoarseNeighborList();
  coarseNeighborList = findNestedNeighbors(mMasterGridCellIndex, mMasterGridLevel);

  // Post conditions.
  ENSURE(coarseNeighborList.size() >= this->masterList().size());
  ENSURE(mMasterGridLevel >= 0 and mMasterGridLevel < mMaxGridLevels);
}

//------------------------------------------------------------------------------
// Set the master list of nodes for a given enter and exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
setMasterList(const GeomPlane<Dimension>& enterPlane,
              const GeomPlane<Dimension>& exitPlane) {

  REQUIRE(valid());
  REQUIRE(enterPlane.parallel(exitPlane));

  // Make references to the two neighbor lists we'll be building for convenience.
  const Field<Dimension, Vector>& nodePositions = this->nodeList().positions();

  // Clear out the master and coarse neighbor lists.
  vector<int>& masterNodeList = this->accessMasterList();
  vector<int>& coarseNeighborList = this->accessCoarseNeighborList();
  masterNodeList = vector<int>();
  coarseNeighborList = vector<int>();

  // Loop over all occupied grid cells.
  int gridLevelID = 0;
  for (typename vector< vector< GridCellIndex<Dimension> > >::const_iterator 
         gridLevelItr = mOccupiedGridCells.begin();
       gridLevelItr < mOccupiedGridCells.end(); ++gridLevelItr, ++gridLevelID) {
    CHECK(gridLevelID >= 0 and gridLevelID < numGridLevels());

    if (mGridLevelOccupied[gridLevelID] == 1) {

      // The size of grid cells on this level.
      const Scalar dxgc = 1.0/mGridCellSizeInv[gridLevelID];

      // Loop over all the grid cells occupied on this grid level.
      for (typename vector< GridCellIndex<Dimension> >::const_iterator
             gridCellItr = (*gridLevelItr).begin();
           gridCellItr < (*gridLevelItr).end(); ++gridCellItr) {

        // Easier to type version of this grid cell.
        const GridCellIndex<Dimension>& gridCell = *gridCellItr;

	// The min and max gridcells on this level that can interact
        // with this grid cell.
	GridCellIndex<Dimension> minGridCell = gridCell - mGridCellInfluenceRadius;
	GridCellIndex<Dimension> maxGridCell = gridCell + mGridCellInfluenceRadius;

        // Find the min and max coordinates in the above grid cell range.
        const Vector rMin = dxgc*minGridCell;
        const Vector rMax = dxgc*(maxGridCell + 1);

        // Min and max positions in the gridcell.
        const Vector rgcMin = dxgc*gridCell;
        const Vector rgcMax = dxgc*(gridCell + 1);

	// If this grid cell is in range of the enter plane it's nodes go into 
	// the control list.
        if ((rMin <= enterPlane and rMax >= enterPlane) or
            (rMin >= enterPlane and rMax <= enterPlane)) {

	  // Loop over the nodes in this grid cell and put them in the control
	  // set.
          int nodeID = headOfGridCell(gridCell, gridLevelID);
          while (nodeID != mEndOfLinkList) {
            CHECK((nodeID >= 0 and nodeID < this->nodeList().numNodes()) or
                  this->nodeList().numNodes() == 0);
            masterNodeList.push_back(nodeID);
            nodeID = nextNodeInCell(nodeID);
            CHECK((nodeID >= 0 and nodeID < this->nodeList().numNodes()) or
                  this->nodeList().numNodes() == 0 or
                  nodeID == mEndOfLinkList);
          }

          // Map the the min and max positions in this grid cell through
          // the enter/exit planes, and find the range of nodes on the exit
          // side that can interact with this grid cell.  They go in the
          // coarse neighbor list.
          Vector r0 = mapPositionThroughPlanes(rgcMin, 
                                               enterPlane, -exitPlane);
          Vector r1 = mapPositionThroughPlanes(rgcMax,
                                               enterPlane, -exitPlane);
          Vector rgcMinExit, rgcMaxExit;
          for (int i = 0; i < Dimension::nDim; ++i) {
            rgcMinExit(i) = std::min(r0(i), r1(i));
            rgcMaxExit(i) = std::max(r0(i), r1(i));
          }

          // Cheat a bit, and push these a little closer together to catch
          // the most typical case where rgcMinExit and rgcMaxExit actually
          // lie on the boundaries of the same gridcell.
          rgcMinExit += 1e-4*dxgc;
          rgcMaxExit -= 1e-4*dxgc;
          CHECK(rgcMinExit <= rgcMaxExit);

          // Determine the grid cells containing these mapped positions.
          const GridCellIndex<Dimension> warpGridCellMin = gridCellIndex(rgcMinExit,
                                                                         gridLevelID);
          const GridCellIndex<Dimension> warpGridCellMax = gridCellIndex(rgcMaxExit,
                                                                         gridLevelID);
          CHECK(warpGridCellMin <= warpGridCellMax);

          // Add the nodes affecting these grid cells (if we haven't already 
          // done so.)
          const vector<int> neighborIndicies = findNestedNeighbors(warpGridCellMin,
                                                                   gridLevelID);
          {
            for (vector<int>::const_iterator idItr = neighborIndicies.begin();
                 idItr < neighborIndicies.end();
                 ++idItr) coarseNeighborList.push_back(*idItr);
          }
          {
            const vector<int> neighborIndicies = findNestedNeighbors(warpGridCellMax,
                                                                     gridLevelID);
            for (vector<int>::const_iterator idItr = neighborIndicies.begin();
                 idItr < neighborIndicies.end();
                 ++idItr) coarseNeighborList.push_back(*idItr);
          }

        } else if ((rMin <= exitPlane and rMax >= exitPlane) or
                   (rMin >= exitPlane and rMax <= exitPlane)) {

	  // If this gridCell overlaps the exit plane, then add its (and its 
          // neighbors) nodes to the coarse neighbor list.
          int nodeID = headOfGridCell(gridCell, gridLevelID);
          while (nodeID != mEndOfLinkList) {
            CHECK((nodeID >= 0 and nodeID < this->nodeList().numNodes()) or
                  this->nodeList().numNodes() == 0);
            coarseNeighborList.push_back(nodeID);
            nodeID = nextNodeInCell(nodeID);
            CHECK((nodeID >= 0 and nodeID < this->nodeList().numNodes()) or
                  this->nodeList().numNodes() == 0 or
                  nodeID == mEndOfLinkList);
          }

          // Map the the min and max positions in this grid cell through
          // the exit/enter planes, and find the range of nodes on the enter
          // side that can interact with this grid cell.  They go in the
          // master neighbor list.
          Vector r0 = mapPositionThroughPlanes(rgcMin, 
                                               exitPlane, -enterPlane);
          Vector r1 = mapPositionThroughPlanes(rgcMax,
                                               exitPlane, -enterPlane);
          Vector rgcMinEnter, rgcMaxEnter;
          for (int i = 0; i < Dimension::nDim; ++i) {
            rgcMinEnter(i) = std::min(r0(i), r1(i));
            rgcMaxEnter(i) = std::max(r0(i), r1(i));
          }

          // Cheat a bit, and push these a little closer together to catch
          // the most typical case where rgcMinExit and rgcMaxExit actually
          // lie on the boundaries of the same gridcell.
          rgcMinEnter += 1e-4*dxgc;
          rgcMaxEnter -= 1e-4*dxgc;
          CHECK(rgcMinEnter <= rgcMaxEnter);

          // Determine the grid cells containing these mapped positions.
          const GridCellIndex<Dimension> warpGridCellMin = gridCellIndex(rgcMinEnter,
                                                                         gridLevelID);
          const GridCellIndex<Dimension> warpGridCellMax = gridCellIndex(rgcMaxEnter,
                                                                         gridLevelID);
          CHECK(warpGridCellMin <= warpGridCellMax);

          // For each grid cell in this range, add it's node to the master 
          // lists.
          {
            const vector<int> neighborIndicies = findNestedNeighbors(warpGridCellMin,
                                                                     gridLevelID);
            for (vector<int>::const_iterator idItr = neighborIndicies.begin();
                 idItr < neighborIndicies.end();
                 ++idItr) masterNodeList.push_back(*idItr);
          }
          {
            const vector<int> neighborIndicies = findNestedNeighbors(warpGridCellMax,
                                                                     gridLevelID);
            for (vector<int>::const_iterator idItr = neighborIndicies.begin();
                 idItr < neighborIndicies.end();
                 ++idItr) masterNodeList.push_back(*idItr);
          }
	}
      }
    }
  }

  // Finally, make sure that we have unique lists of ghost nodes IDs.
  sort(masterNodeList.begin(), masterNodeList.end());
  typename vector<int>::iterator uniqueEnd = unique(masterNodeList.begin(),
                                                    masterNodeList.end());
  masterNodeList.erase(uniqueEnd, masterNodeList.end());
  sort(coarseNeighborList.begin(), coarseNeighborList.end());
  uniqueEnd = unique(coarseNeighborList.begin(),
                     coarseNeighborList.end());
  coarseNeighborList.erase(uniqueEnd, coarseNeighborList.end());

//   // Set the per field coarse data caches for this NodeList.
//   this->nodeList().notifyFieldsCacheCoarseValues();
}

// //------------------------------------------------------------------------------
// // Find any nodes that interact with the given GridCellIndex.
// // This version is based on the oct tree data.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// vector<int>
// NestedGridNeighbor<Dimension>::
// findNestedNeighbors(const GridCellIndex<Dimension>& gridCell, 
// 		    const int gridLevel) const {

//   // Prepare our result.
//   vector<int> result;

//   // The grid cell range for the base grid level.
//   const GC baseMin = gridCell - mGridCellInfluenceRadius;
//   const GC baseMax = gridCell + mGridCellInfluenceRadius;

//   // What grid level should we start with?
//   int currentGridLevel = min(mFirstParentGridLevel, max(0, gridLevel - 1));
//   vector<GC> parents;
//   for (typename OctTree::const_iterator itr = mDaughterCells[currentGridLevel].begin();
//        itr != mDaughterCells[currentGridLevel].end();
//        ++itr) parents.push_back(itr->first);

//   // Iterate down the grid levels until we hit everyone with any contribution.
//   const NeighborSearchType searchType = this->neighborSearchType();
//   GC targetMin, targetMax;
//   while (currentGridLevel < mMaxGridLevels and parents.size() > 0) {
//     vector<GC> newParents;

//     // The range of cells we want on this grid level.
//     translateGridCellRange(baseMin, baseMax, gridLevel, currentGridLevel, targetMin, targetMax);
//     const int radius = (searchType == GatherScatter ? mGridCellInfluenceRadius * intpow2(max(0, currentGridLevel - gridLevel)) :
//                         searchType == Gather ?        mGridCellInfluenceRadius / intpow2(max(0, gridLevel - currentGridLevel)) + 1 :
//                                                       mGridCellInfluenceRadius);
//     targetMin -= radius;
//     targetMax += radius;

//     // Iterate over the current parents, and check if they are in the target region.
//     for (typename vector<GC>::const_iterator parentItr = parents.begin();
//          parentItr != parents.end();
//          ++parentItr) {
//       if (parentItr->inRange(targetMin, targetMax)) {
//         appendNodesInCell(*parentItr, currentGridLevel, result);
//         const vector<GC>& daughters = daughterCells(*parentItr, currentGridLevel);
//         copy(daughters.begin(), daughters.end(), back_inserter(newParents));
//       }
//     }
 
//     // The next generation.
//     parents = newParents;
//     ++currentGridLevel;
//   }

//   // That's it.
//   return result;
// }

//------------------------------------------------------------------------------
// Find any nodes that interact with the given GridCellIndex.
// This version just searches the occupied grid cells level by level.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
NestedGridNeighbor<Dimension>::
findNestedNeighbors(const GridCellIndex<Dimension>& gridCell, 
		    int gridLevel) const {

  // Vector of node indicies which we're going to build up.
  vector<int> result;

  // The grid cell range for the base grid level.
  const GC baseMin = gridCell - mGridCellInfluenceRadius;
  const GC baseMax = gridCell + mGridCellInfluenceRadius;

  // Loop over all occupied grid levels.
  const NeighborSearchType searchType = this->neighborSearchType();
  GC targetMin, targetMax;
  for (int currentGridLevel = 0; 
       currentGridLevel != mMaxGridLevels;
       ++currentGridLevel) {
    if (mGridLevelOccupied[currentGridLevel] == true) {

      // Find the range of grid cells on this grid level that overlap
      // this gridcells potential neighbor influence.
      translateGridCellRange(baseMin, baseMax, gridLevel, currentGridLevel, targetMin, targetMax);
      const int radius = (searchType == GatherScatter ? mGridCellInfluenceRadius * intpow2(max(0, currentGridLevel - gridLevel)) :
			  searchType == Gather ?        mGridCellInfluenceRadius / intpow2(max(0, gridLevel - currentGridLevel)) + 1 :
                                                        mGridCellInfluenceRadius);
      targetMin -= radius;
      targetMax += radius;
      CHECK(targetMin <= targetMax);

      // Loop over the occupied cells in the neighbor grid cell range,
      // and add their nodes to the result.
      vector<GC> gridCells;
      occupiedGridCellsInRange(gridCells, targetMin, targetMax, currentGridLevel);
      for (typename vector<GC>::iterator gcItr = gridCells.begin();
	   gcItr != gridCells.end();
	   ++gcItr) appendNodesInCell(*gcItr, currentGridLevel, result);
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Set the number of grid levels.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
numGridLevels(const int numGridLevels) {
  VERIFY(numGridLevels > 0 and numGridLevels < 32);
  mMaxGridLevels = numGridLevels;

  // Resize the occupied grid level array.
  mGridLevelOccupied.resize(numGridLevels);

  // If there is enough information, recalculate all the grid info.
  double topGridCellSize = topGridSize();
  if (topGridCellSize > 0.0 and this->kernelExtent() > 0.0) {
      mGridLevelConst0 = log(topGridCellSize*mGridCellInfluenceRadius)*ln2inverse;

      // Build the vector of gridcell sizes.
      mGridCellSizeInv.resize(mMaxGridLevels);
      for (unsigned i = 0; i != numGridLevels; ++i) {
        mGridCellSizeInv[i] = double(intpow2(i))/topGridCellSize;
      }

      // Rebuild the node per gridcell info.
      updateNodes();

      ENSURE(valid());
  }
}

//------------------------------------------------------------------------------
// Set the top grid cell size.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::topGridSize(const double topGridSize) {
  REQUIRE(topGridSize > 0.0);
  REQUIRE(mMaxGridLevels > 0);

  // Build the vector of gridcell sizes.
  for (unsigned i = 0; i != mMaxGridLevels; ++i) {
    mGridCellSizeInv[i] = double(intpow2(i))/topGridSize;
  }

  // We have to recalculate the grid level constant.
  mGridLevelConst0 = log(topGridSize*mGridCellInfluenceRadius)*ln2inverse;

  // If we have enough info, assign nodes to the grid cells.
  if (readyToAssignNodes()) updateNodes();
}

//------------------------------------------------------------------------------
// Return the list of occupied grid cells in the given range.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NestedGridNeighbor<Dimension>::
// occupiedGridCellsInRange(vector< GridCellIndex<Dimension> >& gridCells,
//                          const GridCellIndex<Dimension>& minGridCell,
//                          const GridCellIndex<Dimension>& maxGridCell,
//                          const int gridLevel) const {

//   REQUIRE(minGridCell <= maxGridCell);
//   REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);

//   // Clear any existing result.
//   gridCells = vector<GC>();

//   // What grid level should we start with?
//   int currentGridLevel = min(mFirstParentGridLevel, max(0, gridLevel - 1));
//   GC targetMin, targetMax;
//   translateGridCellRange(minGridCell, maxGridCell, gridLevel, currentGridLevel, targetMin, targetMax);
//   for (typename OctTree::const_iterator itr = mDaughterCells[currentGridLevel].begin();
//        itr != mDaughterCells[currentGridLevel].end();
//        ++itr) {
//     if (itr->first.inRange(targetMin, targetMax)) gridCells.push_back(itr->first);
//   }

//   // Iterate down the grid levels until we hit the desired level.
//   while (currentGridLevel < gridLevel and gridCells.size() > 0) {
//     ++currentGridLevel;
//     vector<GC> newParents;

//     // The range of cells we want on this grid level.
//     translateGridCellRange(minGridCell, maxGridCell, gridLevel, currentGridLevel, targetMin, targetMax);

//     // Iterate over the current parents, and check if they are in the target region.
//     for (typename vector<GC>::const_iterator parentItr = gridCells.begin();
//          parentItr != gridCells.end();
//          ++parentItr) {
//       const vector<GC>& daughters = daughterCells(*parentItr, currentGridLevel);
//       for (typename vector<GC>::const_iterator itr = daughters.begin();
//            itr != daughters.end();
//            ++itr) {
//         if (itr->inRange(targetMin, targetMax)) newParents.push_back(*itr);
//       }
//     }

//     // The next generation.
//     gridCells = newParents;
//   }
// }

template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
occupiedGridCellsInRange(vector< GridCellIndex<Dimension> >& gridCells,
                         const GridCellIndex<Dimension>& minGridCell,
                         const GridCellIndex<Dimension>& maxGridCell,
                         const int gridLevelID) const {

  REQUIRE(minGridCell <= maxGridCell);

  // Count how many grid cells are in the given range.
  const int numCellsInRange = (maxGridCell - minGridCell + 1).productElements();

  // Determine whether it's less work to check every grid cell in the range,
  // or to scan every occupied cell on this level.
  if (occupiedGridCells(gridLevelID).size() < numCellsInRange) {

    // Loop over the occupied grid cells on the given grid level.
    typename vector<GC>::const_iterator begin = occupiedGridCells(gridLevelID).begin();
    typename vector<GC>::const_iterator end = occupiedGridCells(gridLevelID).end();
    gridCells.reserve(occupiedGridCells(gridLevelID).size());
    for (typename vector<GC>::const_iterator gridCellItr = begin; gridCellItr != end; ++gridCellItr) {

      // If this grid cell is in the required range, add it to the list.
      if ((*gridCellItr).inRange(minGridCell, maxGridCell)) gridCells.push_back(*gridCellItr);
    }

  } else {

    // Scan all the grid cells in the given range.
    gridCells.reserve(numCellsInRange);
    for (GC gridCell = minGridCell; gridCell <= maxGridCell; incrementGridCell(gridCell, minGridCell, maxGridCell)) {
      const int nodeID = headOfGridCell(gridCell, gridLevelID);
      if (nodeID != mEndOfLinkList) gridCells.push_back(gridCell);
    }
  }
}

//------------------------------------------------------------------------------
// Map a grid cell through an enter/exit plane combo.
// At the moment this only works correctly for the case where we are mapping
// along one of the cardinal directions.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellIndex<Dimension>
NestedGridNeighbor<Dimension>::
mapGridCell(const GridCellIndex<Dimension>& gridCell,
	    const int gridLevel,
            const GeomPlane<Dimension>& enterPlane0,
            const GeomPlane<Dimension>& exitPlane0) const {
  REQUIRE(gridLevel >= 0 and gridLevel < mMaxGridLevels);
  REQUIRE(enterPlane0.parallel(exitPlane0));

  // Map the gridCell to a Vector position, and likewise convert the points of
  // the planes to the appropriate coordinates to work on this gridlevel.
  Vector r;
  for (int i = 0; i < Dimension::nDim; ++i) r(i) = gridCell(i);
  Vector renter = enterPlane0.point()*mGridCellSizeInv[gridLevel];
  Vector rexit = exitPlane0.point()*mGridCellSizeInv[gridLevel];
  const GeomPlane<Dimension> enterPlane(renter, enterPlane0.normal());
  const GeomPlane<Dimension> exitPlane(rexit, exitPlane0.normal());

  // Now map the continuous version of the gridcell just like we do in 
  // PlanarBoundary::mapPosition.
  const Vector deltaEnter = (r - enterPlane.point()).dot(enterPlane.normal())*enterPlane.normal();
  const Vector deltaExit = (r - exitPlane.point()).dot(exitPlane.normal())*exitPlane.normal();
  const Vector rmap = r - deltaExit + deltaEnter.magnitude()*exitPlane.normal();

  // Now convert this continuous mapped postion back to a discrete GridCellIndex.
  GridCellIndex<Dimension> result;
  for (int i = 0; i < Dimension::nDim; ++i) result(i) = int(rmap(i));

//   cerr << gridCell << " "
//        << deltaExit << " "
//        << deltaEnter << " "
//        << exitPlane << " "
//        << rmap << " "
//        << result
//        << endl;
  return result;
}

//------------------------------------------------------------------------------
// Determine if the NestedGridNeighbor is in a minimally valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NestedGridNeighbor<Dimension>::valid() const {
  const NodeList<Dimension>& nodes = this->nodeList();

  if (mMaxGridLevels <= 0) {
    cerr << "NestedGridNeighbor::valid: invalid mMaxGridLevels "
         << mMaxGridLevels << endl;
    return false;
  }

  if (topGridSize() <= 0.0) {
    cerr << "NestedGridNeighbor::valid: invalid topGridSize "
         << topGridSize() << endl;
    return false;
  }

  if (this->kernelExtent() <= 0.0) {
    cerr << "NestedGridNeighbor::valid: invalid kernelExtent "
         << this->kernelExtent() << endl;
    return false;
  }

  if (!fuzzyEqual(mGridLevelConst0, log(topGridSize()*mGridCellInfluenceRadius)*ln2inverse)) {
    cerr << "NestedGridNeighbor::valid: invalid mGridLevelConst0 : "
                   << mGridLevelConst0 << " "
                   << log(topGridSize()*mGridCellInfluenceRadius)*ln2inverse << endl;
    return false;
  }

  if (mGridCellSizeInv.size() != mMaxGridLevels) {
    cerr << "NestedGridNeighbor::valid: invalid mMaxGridLevels  "
         << mMaxGridLevels << endl;
    return false;
  }

  if (mGridLevelOccupied.size() != mMaxGridLevels) {
    cerr << "NestedGridNeighbor::valid: invalid mGridLevelOccupied "
         << mGridLevelOccupied.size() << endl;
    return false;
  }

  // Check that each node is assigned to one and only one grid cell.
  {
    Field<Dimension, int> count("count", nodes);
    for (typename vector<SafeIndexMap<GridCellIndex<Dimension>, int> >::const_iterator itr1 = mGridCellHead.begin();
         itr1 != mGridCellHead.end();
         ++itr1) {
      for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr2 = itr1->begin();
           itr2 != itr1->end();
           ++itr2) {
        int i = itr2->second;
        while (i != mEndOfLinkList) {
          count(i)++;
          i = mNextNodeInCell[i];
        }
      }
    }
    for (int i = 0; i != nodes.numInternalNodes(); ++i) {
      if (count(i) != 1) {
        cerr << "NestedGridNeighbor::valid : incorrect count of assignment to gridcell for node "
             << i << " " 
             << count(i) << endl;
        return false;
      }
    }
  }

//   // Check that each node is reachable through the oct tree.
//   {
//     vector<int> nodeIDs;

//     // Initialize from the top level.
//     vector<GC> parents;
//     for (typename OctTree::const_iterator parentItr = mDaughterCells[0].begin();
//          parentItr != mDaughterCells[0].end();
//          ++parentItr) {
//       const GC& parent = parentItr->first;
//       appendNodesInCell(parent, 0, nodeIDs);
//       const vector<GC>& daughters = parentItr->second;
//       CHECK(daughters.size() <= intpow2(Dimension::nDim));
//       copy(daughters.begin(), daughters.end(), back_inserter(parents));
//     }
      
//     // Descend through the grid levels until we're done.
//     int gridLevel = 1;
//     while (gridLevel < mMaxGridLevels and parents.size() > 0) {
//       vector<GC> newParents;
//       for (typename vector<GC>::const_iterator parentItr = parents.begin();
//            parentItr != parents.end();
//            ++parentItr) {
//         appendNodesInCell(*parentItr, gridLevel, nodeIDs);
//         const vector<GC>& daughters = daughterCells(*parentItr, gridLevel);
//         copy(daughters.begin(), daughters.end(), back_inserter(newParents));
//       }
//       ++gridLevel;
//       parents = newParents;
//     }

//     // Now we can check if all nodes are accounted for.
//     sort(nodeIDs.begin(), nodeIDs.end());
//     if (nodeIDs.size() != nodes.numNodes()) {
//       cerr << "NestedGridNeighbor::valid : oct tree does not hit the correct number of nodes : " 
//            << nodeIDs.size() << " "
//            << nodes.numNodes()  << " "
//            << endl;
//       return false;
//     }
//     for (int i = 0; i != nodes.numNodes(); ++i) {
//       if (nodeIDs[i] != i) {
//         cerr << "NestedGridNeighbor::valid : oct tree node indicies wrong?" << endl;
//         return false;
//       }
//     }
//   }

//   // Check that each node can be found down the correct branch of the tree.
//   {
//     for (int i = 0; i != nodes.numNodes(); ++i) {
//       const int gridLeveli = gridLevel(i);

//       for (int gridLevel = 0; gridLevel != gridLeveli; ++gridLevel) {
//         const GC gcParent = gridCellIndex(i, gridLevel);
//         const GC gcDaughter = gridCellIndex(i, gridLevel + 1);
//         const vector<GC>& daughters = daughterCells(gcParent, gridLevel);
//         if (not (cellInTree(gcParent, gridLevel) and 
//                  ((gridLevel + 1 == gridLeveli) or cellInTree(gcDaughter, gridLevel + 1)) and
//                  (find(daughters.begin(), daughters.end(), gcDaughter) != daughters.end()))) {
//           cerr << "Path down oct tree to node " << i << " is broken." 
//                << endl
//                << "Node in grid cell " << gridCellIndex(i, gridLeveli)
//                << " on grid level " << gridLeveli 
//                << endl;
//           cerr << "Test conditions : "
//                << cellInTree(gcParent, gridLevel) << " "
//                << ((gridLevel + 1 == gridLeveli) or cellInTree(gcDaughter, gridLevel + 1)) << " "
//                << (find(daughters.begin(), daughters.end(), gcDaughter) != daughters.end()) 
//                << endl
//                << " Parent tree thus far: "
//                << endl;
//           for (int j = 0; j != gridLevel; ++j) {
//             cerr << "    "
//                  << j << " "
//                  << gridCellIndex(i, j) << endl;
//           }
//           cerr << "Current grid level : " << gridLevel
//                << "          gcParent : " << gcParent 
//                << "        gcDaughter : " << gcDaughter
//                << endl
//                << "Full set of daughters : "
//                << endl;
//           for (typename vector<GC>::const_iterator itr = daughters.begin();
//                itr != daughters.end();
//                ++itr) cerr << "    " << *itr << endl;
//           return false;
//         }
//       }
//     }
//   }

  return true;
}

//------------------------------------------------------------------------------
// Check whether we have enough info to assign nodes to grid cells.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NestedGridNeighbor<Dimension>::readyToAssignNodes() const {
  if (mMaxGridLevels == 0) return false;
  return (this->kernelExtent() > 0.0 and
          fuzzyEqual(mGridLevelConst0, log(topGridSize()*mGridCellInfluenceRadius)*ln2inverse) and
          topGridSize() > 0.0 and
          mGridLevelOccupied.size() == numGridLevels());
}

//------------------------------------------------------------------------------
// Build the node per grid info from scratch.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::updateNodes() {
  REQUIRE(readyToAssignNodes());

  // Resize things to the correct sizes.
  const int numNodes = this->nodeList().numNodes();
  for (int gridLevelID = 0; gridLevelID < numGridLevels(); ++gridLevelID) {
    CHECK(gridLevelID < mGridLevelOccupied.size());
    CHECK(gridLevelID < mGridCellHead.size());
    CHECK(gridLevelID < mNodeInCell.size());
    mGridLevelOccupied[gridLevelID] = false;
    mGridCellHead[gridLevelID] = SafeIndexMap<GridCellIndex<Dimension>, int>();
    mNodeInCell[gridLevelID].resize(numNodes);
    mOccupiedGridCells[gridLevelID] = vector< GridCellIndex<Dimension> >();
  }
  mNextNodeInCell.resize(numNodes);
  mNodeOnGridLevel.resize(numNodes);

  // Blank out the link list of nodes.
  for (int nodeID = 0; nodeID != numNodes; ++nodeID) {
    CHECK(nodeID < mNextNodeInCell.size());
    mNextNodeInCell[nodeID] = mEndOfLinkList;
  }

  // Loop over all the nodes and assign them their home grid level & cell.
  for (int nodeID = 0; nodeID != numNodes; ++nodeID) {

    // Determine the gridlevel and gridcell for this node.
    const int gridLevelID = gridLevel(nodeID);
    const GridCellIndex<Dimension> gridCellID = gridCellIndex(nodeID, gridLevelID);
    CHECK(gridLevelID >= 0 and gridLevelID < mMaxGridLevels);

    // Assign the nodes grid level.
    CHECK(nodeID >= 0 and nodeID < mNodeOnGridLevel.size());
    mNodeOnGridLevel[nodeID] = gridLevelID;

    // Link this node into the grid cell/grid level.
    linkNode(nodeID, gridLevelID, gridCellID);

    // Flag the grid level this node is on as occupied.
    CHECK(gridLevelID >= 0 and gridLevelID < mGridLevelOccupied.size());
    mGridLevelOccupied[gridLevelID] = true;
  }

  // Now we need to identify which cell on each grid level each node occupies.
  for (int gridLevelID = 0; gridLevelID != mMaxGridLevels; ++gridLevelID) {
    for (int nodeID = 0; nodeID != numNodes; ++nodeID) {
      CHECK(gridLevelID >= 0 and gridLevelID < mNodeInCell.size());
      CHECK(nodeID >= 0 and nodeID < mNodeInCell[gridLevelID].size());
      mNodeInCell[gridLevelID][nodeID] = gridCellIndex(nodeID, gridLevelID);
    }
  }

  // Create the list of occupied grid cells.
  rebuildOccupiedGridCells();

//   // Build the tree.
//   rebuildOctTree();

  // Force the node extents to be calculated.
  this->setNodeExtents();
  for (int nodeID = 0; nodeID != numNodes; ++nodeID) {
    CHECK(fuzzyLessThanOrEqual(this->nodeExtentField()(nodeID).maxElement(),
                                topGridSize()));
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    // Make sure each node is only linked into one grid cell.
    const NodeList<Dimension>& nodes = this->nodeList();
    Field<Dimension, uint64_t> count("count check", nodes);
    const vector<int> occupiedGLs = occupiedGridLevels();
    for (vector<int>::const_iterator glItr = occupiedGLs.begin();
         glItr != occupiedGLs.end();
         ++glItr) {
      const vector<GridCellIndex<Dimension> >& occupiedCells = occupiedGridCells(*glItr);
      for (typename vector<GridCellIndex<Dimension> >::const_iterator gcItr = occupiedCells.begin();
           gcItr != occupiedCells.end();
           ++gcItr) {
        int i = headOfGridCell(*gcItr, *glItr);
        while (i != mEndOfLinkList) {
          count(i) += 1;
          i = nextNodeInCell(i);
        }
      }
    }
    for (int i = 0; i != nodes.numNodes(); ++i) ENSURE(count(i) == 1);
  }
  END_CONTRACT_SCOPE;

}

//------------------------------------------------------------------------------
// Reassign the given nodes to gridcells.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
updateNodes(const vector<int>& nodeIDs) {

  REQUIRE(readyToAssignNodes());

  // Make sure the internal data structures are sized correctly.
  int numNodes = this->nodeList().numNodes();
  for (int gridLevelID = 0; gridLevelID < numGridLevels(); ++gridLevelID) {
    REQUIRE(gridLevelID < mNodeInCell.size());
    mNodeInCell[gridLevelID].resize(numNodes);
  }
  mNextNodeInCell.resize(numNodes);
  mNodeOnGridLevel.resize(numNodes);

  // Loop over the nodes we're working on.
  for (typename vector<int>::const_iterator nodeIDItr = nodeIDs.begin();
       nodeIDItr < nodeIDs.end(); ++nodeIDItr) {
    CHECK(*nodeIDItr >= 0 and *nodeIDItr < this->nodeList().numNodes());

    // Determine the nodes current grid level and cell.
    const int gridLevelID = gridLevel(*nodeIDItr);
    CHECK(gridLevelID >= 0 and gridLevelID < mMaxGridLevels);
    const GridCellIndex<Dimension> gridCellID = gridCellIndex(*nodeIDItr, gridLevelID);

    // If the new grid cell doesn't match the old, unlink and reassign.
//     if (gridLevelID != mNodeOnGridLevel[*nodeIDItr] or
//         gridCellID != mNodeInCell[mNodeOnGridLevel[*nodeIDItr]][*nodeIDItr]) {

      // Unlink this node from it's original grid cell.
      unlinkNode(*nodeIDItr, 
                 mNodeOnGridLevel[*nodeIDItr],
                 mNodeInCell[mNodeOnGridLevel[*nodeIDItr]][*nodeIDItr]);

      // Link the node into it's new grid cell.
      linkNode(*nodeIDItr, gridLevelID, gridCellID);

      // Assign the node's new grid level.
      mNodeOnGridLevel[*nodeIDItr] = gridLevelID;
//     }

    // For all gridlevels we need to identify the grid cell that contains this
    // node.
    for (int gridLevelID = 0; gridLevelID < mMaxGridLevels; ++gridLevelID) {
      CHECK(gridLevelID >= 0 and gridLevelID < mNodeInCell.size());
      CHECK(*nodeIDItr >= 0 and *nodeIDItr < mNodeInCell[gridLevelID].size());
      mNodeInCell[gridLevelID][*nodeIDItr] = gridCellIndex(*nodeIDItr, gridLevelID);
    }

    // Flag the grid level this node is on as occupied.
    CHECK(gridLevelID >= 0 and gridLevelID < mGridLevelOccupied.size());
    mGridLevelOccupied[gridLevelID] = true;
  }

  // Create the list of occupied grid cells.
  rebuildOccupiedGridCells();

//   // Build the tree.
//   rebuildOctTree();

  // Force the node extents to be calculated.
  this->setNodeExtents(nodeIDs);
  for (int nodeID = 0; nodeID < this->nodeList().numNodes(); ++nodeID) {
    CHECK(fuzzyLessThanOrEqual(this->nodeExtentField()(nodeID).maxElement(),
                                topGridSize()));
  }

}

//------------------------------------------------------------------------------
// Build the list of occupied grid cells from scratch.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridNeighbor<Dimension>::
rebuildOccupiedGridCells() {

  REQUIRE(readyToAssignNodes());
  REQUIRE(mGridCellHead.size() == numGridLevels());
  REQUIRE(mOccupiedGridCells.size() == numGridLevels());

  // Loop over all grid levels.
  for (int gridLevelID = 0; gridLevelID < numGridLevels(); ++gridLevelID) {

    // Erase any old info.
    mOccupiedGridCells[gridLevelID] = vector< GridCellIndex<Dimension> >();

    // Loop over all the grid cell heads for this grid level.
    if (!(mGridCellHead[gridLevelID].empty())) {
      for (typename map<GridCellIndex<Dimension>, int>::const_iterator
             headItr = mGridCellHead[gridLevelID].begin();
           headItr != mGridCellHead[gridLevelID].end(); ++headItr) {

        // Add this grid cell to the list for this grid level.
        CHECK((*headItr).second >= 0 and (*headItr).second < this->nodeList().numNodes());
        mOccupiedGridCells[gridLevelID].push_back((*headItr).first);
      }
    }
    CHECK(mOccupiedGridCells[gridLevelID].size() == mGridCellHead[gridLevelID].size());
  }

}

//------------------------------------------------------------------------------
// Build the oct tree from scratch.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NestedGridNeighbor<Dimension>::
// rebuildOctTree() {

//   // Clear any pre-exising info.
//   mDaughterCells = vector<OctTree>(mMaxGridLevels);
  
//   // Iterate over each node and assign it through the tree down to its home
//   // grid level.
//   const int numNodes = this->nodeList().numNodes();
//   for (int i = 0; i != numNodes; ++i) {
//     const int gridLeveli = gridLevel(i);

//     // Iterate over the grid levels until we get to the home level for this
//     // node.
//     for (int gridLevel = 0; gridLevel != gridLeveli; ++gridLevel) {
//       const GC parent = gridCellIndex(i, gridLevel);
//       const GC daughter = gridCellIndex(i, gridLevel + 1);
//       vector<GC>& siblings = mDaughterCells[gridLevel].unsafeIndex(parent);
//       siblings.reserve(intpow2(Dimension::nDim));
//       typename vector<GC>::iterator itr = lower_bound(siblings.begin(), siblings.end(), daughter);
//       if (itr == siblings.end() or *itr != daughter) siblings.insert(itr, daughter);
//       CHECK( lower_bound(siblings.begin(), siblings.end(), daughter) != siblings.end() and
//             *lower_bound(siblings.begin(), siblings.end(), daughter) == daughter);
//     }
//   }

//   // Find the last grid level with only one cell acting as a parent.  This will
//   // be our starting point for walking the tree.
//   mFirstParentGridLevel = 0;
//   while (mFirstParentGridLevel < mMaxGridLevels and 
//          mDaughterCells[mFirstParentGridLevel].size() > 0 and
//          mDaughterCells[mFirstParentGridLevel].size() < 2 and
//          mOccupiedGridCells[mFirstParentGridLevel].size() == 0) ++mFirstParentGridLevel;
//   mFirstParentGridLevel = max(0, mFirstParentGridLevel - 1);
//   ENSURE(mFirstParentGridLevel >= 0 and mFirstParentGridLevel < mMaxGridLevels);

// }

//------------------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------------------
template<typename Dimension>
const double NestedGridNeighbor<Dimension>::ln2inverse = 1.0/log(2.0);
template<typename Dimension>
const int NestedGridNeighbor<Dimension>::mEndOfLinkList = -1;
template<typename Dimension>
const int NestedGridNeighbor<Dimension>::mGridNormalMagnitude = 1024;
}
}

