//---------------------------------Spheral++----------------------------------//
// TreeNeighbor
//
// An SPH neighbor finder based on an oct-tree like structure with specialized
// cell membership criteria.
//
// Based on the algorithm originally described in 
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by J. Michael Owen, Fri Mar 23 15:50:35 PDT 2012
//----------------------------------------------------------------------------//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "TreeNeighbor.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Geometry/Dimension.hh"
#include "DBC.hh"

namespace Spheral {
namespace NeighborSpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeNeighbor<Dimension>::
TreeNeighbor(NodeList<Dimension>& nodeList,
             const NeighborSearchType searchType,
             const double kernelExtent):
  Neighbor<Dimension>(nodeList, searchType, kernelExtent),
  mBoxLength(0.0),
  mGridLevelConst0(0.0),
  mTree() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeNeighbor<Dimension>::
~TreeNeighbor() {
}

//------------------------------------------------------------------------------
// Set the master list using either a scalar or tensor H.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
setMasterList(const Vector& position,
              const Scalar& H) {
  REQUIRE(H > 0.0);
  const Scalar h = 1.0/H;
  this->setTreeMasterList(position, h);
}

template<typename Dimension>
void
TreeNeighbor<Dimension>::
setMasterList(const Vector& position,
              const SymTensor& H) {
  REQUIRE(H.Determinant() > 0.0);
  const Vector hinvValues = H.eigenValues();
  CHECK(hinvValues.minElement() > 0.0);
  const Scalar h = 1.0/hinvValues.minElement();
  this->setTreeMasterList(position, h);
}

template<typename Dimension>
void
TreeNeighbor<Dimension>::
setMasterList(const Vector& position) {
  this->setTreeMasterList(position, 1.0e-30*mBoxLength);
}

//------------------------------------------------------------------------------
// Set the refine neighbor list using either a scalar or tensor H.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
setRefineNeighborList(const Vector& position,
                      const Scalar& H) {
  REQUIRE(H > 0.0);
  this->setTreeRefineNeighborList(position, H*SymTensor::one);
}

template<typename Dimension>
void
TreeNeighbor<Dimension>::
setRefineNeighborList(const Vector& position,
                      const SymTensor& H) {
  REQUIRE(H.Determinant() > 0.0);
  this->setTreeRefineNeighborList(position, H);
}

template<typename Dimension>
void
TreeNeighbor<Dimension>::
setRefineNeighborList(const Vector& position) {
  this->setTreeRefineNeighborList(position, 1.0e30*mBoxLength*SymTensor::one);
}

//------------------------------------------------------------------------------
// Set the master list based on proximity to planes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
setMasterList(const GeomPlane<Dimension>& enterPlane,
              const GeomPlane<Dimension>& exitPlane) {
  // FIXME
}

//------------------------------------------------------------------------------
// Update the internal data for all nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
updateNodes() {

  // Clear our internal data.
  mTree = Tree();

  // Grab the NodeList state.
  const NodeList<Dimension>& nodes = this->nodeList();
  const Field<Dimension, Vector>& positions = nodes.positions();
  const Field<Dimension, SymTensor>& H = nodes.Hfield();

  // Recompute the current box size.
  Vector xmin, xmax;
  globalBoundingBox(positions, xmin, xmax, false);
  mBoxLength = (xmax - xmin).maxElement();
  CHECK(mBoxLength > 0.0);
  mGridLevelConst0 = log(mBoxLength)/log(2.0);

  // Walk all the internal nodes and add them to the tree.
  const size_t n = nodes.numInternalNodes();
  for (unsigned i = 0; i != n; ++i) {
    this->addNodeToTree(positions(i), H(i), i);
  }

  // Set the daughter pointers.
  constructDaughterPtrs(mTree);
}

//------------------------------------------------------------------------------
// Update the internal data for a subset of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
updateNodes(const vector<int>& ndoeIDs) {
  this->updateNodes();  // Punt and just rebuild everything.
}

//------------------------------------------------------------------------------
// Find the grid level for the given smoothing scale.
//------------------------------------------------------------------------------
// Argument units of length
template<typename Dimension>
unsigned
TreeNeighbor<Dimension>::
gridLevel(const double& h) const {   
  const unsigned result = std::max(0U, 
                                   std::min(num1dbits - 1U,
                                            unsigned(mGridLevelConst0 - log(h)/log(2.0))));
  ENSURE(fuzzyLessThanOrEqual(mBoxLength/(1U << result), h*this->kernelExtent(), 1.0e-10) and
         fuzzyGreaterThanOrEqual(mBoxLength/(1U << (result + 1U)), h*this->kernelExtent(), 1.0e-10));
  return result;
}

// Argument units of 1/length
template<typename Dimension>
unsigned
TreeNeighbor<Dimension>::
gridLevel(const typename Dimension::SymTensor& H) const {   
  REQUIRE(H.Determinant() > 0.0);
  return this->gridLevel(1.0/H.eigenValues().minElement());
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTree(const bool globalTree) const {
  stringstream ss;
  CellKey key, ix, iy, iz;
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
  unsigned nlevels = mTree.size();
  if (globalTree) nlevels = allReduce(nlevels, MPI_MAX, MPI_COMM_WORLD);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mTree[ilevel].size());
    if (ilevel < mTree.size()) {
      for (typename TreeLevel::const_iterator itr = mTree[ilevel].begin();
           itr != mTree[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef USE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      extractCellIndices(cell.key, ix, iy, iz);
      ss << "    Cell key=" << cell.key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << ")\n"
         << "         daughters = ( ";
      for (typename vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) ss << *ditr << " ";
      ss << ")\n"
         << "         nodes = [";
      for (unsigned k = 0; k != cell.members.size(); ++k) ss << cell.members[k] << " ";
      ss <<" ]\n";
    }
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeNeighbor<Dimension>::
dumpTreeStatistics(const bool globalTree) const {
  stringstream ss;
  CellKey key, ix, iy, iz;
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
  unsigned nlevels = mTree.size();
  if (globalTree) nlevels = allReduce(nlevels, MPI_MAX, MPI_COMM_WORLD);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mTree[ilevel].size());
    if (ilevel < mTree.size()) {
      for (typename TreeLevel::const_iterator itr = mTree[ilevel].begin();
           itr != mTree[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef USE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    // Count up how much of everything we have.
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    unsigned nparticles = 0;
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      nparticles += cell.members.size();
    }
    ss << "         : nparts = " << nparticles << "\n";
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// Serialize to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
serialize(std::vector<char>& buffer) const {
  packElement(mBoxLength, buffer);
  packElement(mGridLevelConst0, buffer);
  const unsigned nlevels = mTree.size();
  packElement(nlevels, buffer);
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    const unsigned ncells = mTree[ilevel].size();
    packElement(ncells, buffer);
    for (typename TreeLevel::const_iterator itr = mTree[ilevel].begin();
         itr != mTree[ilevel].end();
         ++itr) {
      packElement(itr->first, buffer);
      serialize(itr->second, buffer);
    }
  }
}

//------------------------------------------------------------------------------
// Serialize a cell to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
serialize(const TreeNeighbor<Dimension>::Cell& cell,
          std::vector<char>& buffer) const {
  packElement(cell.key, buffer);
  packElement(cell.daughters, buffer);
  packElement(cell.members, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a tree from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
deserialize(vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) {
  unpackElement(mBoxLength, bufItr, endItr);
  unpackElement(mGridLevelConst0, bufItr, endItr);
  unsigned nlevels, ncells;
  CellKey key;
  Cell cell;
  unpackElement(nlevels, bufItr, endItr);
  mTree.resize(nlevels);
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    unsigned ncells;
    unpackElement(ncells, bufItr, endItr);
    for (unsigned i = 0; i != ncells; ++i) {
      cell.daughters = vector<CellKey>();
      cell.members = vector<int>();
      unpackElement(key, bufItr, endItr);
      deserialize(cell, bufItr, endItr);
      mTree[ilevel][key] = cell;
    }
  }
  constructDaughterPtrs(mTree);
}

//------------------------------------------------------------------------------
// Deserialize a cell from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
deserialize(typename TreeNeighbor<Dimension>::Cell& cell,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unpackElement(cell.key, bufItr, endItr);
  unpackElement(cell.daughters, bufItr, endItr);
  unpackElement(cell.members, bufItr, endItr);
}

//------------------------------------------------------------------------------
// Set the master & coarse neighbor sets by walking the tree.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
setTreeMasterList(const typename Dimension::Vector& position,
                  const double& h) {

  // TAU timers.
  TAU_PROFILE("TreeNeighbor", "::setTreeMasterList", TAU_USER);

  // Set the working master grid level and cell.
  CellKey masterKey, ix_master, iy_master, iz_master;
  const LevelKey masterLevel = gridLevel(h);
  buildCellKey(masterLevel, position, masterKey, ix_master, iy_master, iz_master);
  CHECK(masterLevel >= 0 and masterLevel < num1dbits);

  // Set the master list.
  vector<int>& masterList = this->accessMasterList();
  typename TreeLevel::const_iterator masterItr = mTree[masterLevel].find(masterKey);
  CHECK(masterItr != mTree[masterLevel].end());
  masterList = masterItr->second.members;

  // Find all the potential neighbors.
  vector<int>& coarseNeighborList = this->accessCoarseNeighborList();
  coarseNeighborList = this->findTreeNeighbors(masterLevel, ix_master, iy_master, iz_master);

  // Post conditions.
  ENSURE(coarseNeighborList.size() >= this->masterList().size());
}

//------------------------------------------------------------------------------
// Produce the refined list of potential neighbors for a single node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TreeNeighbor<Dimension>::
setTreeRefineNeighborList(const typename Dimension::Vector& position,
                          const typename Dimension::SymTensor& H) {

  // TAU timers.
  TAU_PROFILE("TreeNeighbor", "::setTreeRefineNeighborList", TAU_USER);

  // Determine the maximum extent of this H tensor in each dimension.
  const Vector extent = this->HExtent(H, this->kernelExtent());
  const Vector minExtent = position - extent;
  const Vector maxExtent = position + extent;

  // Use precull to set the refined neighbor list.
  const std::vector<int>& coarseList = this->coarseNeighborList();
  std::vector<int>& refineList = this->accessRefineNeighborList();
  refineList = this->precullList(position, position, minExtent, maxExtent, coarseList);
}

//------------------------------------------------------------------------------
// Walk the tree to find the neighbors for the given (gridlevel, gridcell) 
// pair.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
TreeNeighbor<Dimension>::
findTreeNeighbors(const LevelKey& masterLevel,
                  const typename TreeNeighbor<Dimension>::CellKey& ix_master,
                  const typename TreeNeighbor<Dimension>::CellKey& iy_master,
                  const typename TreeNeighbor<Dimension>::CellKey& iz_master) const {

  // TAU timers.
  TAU_PROFILE("TreeNeighbor", "::findTreeNeighbors", TAU_USER);

  // Declare variables.
  LevelKey ilevel = 0;
  CellKey ix_min, iy_min, iz_min, ix_max, iy_max, iz_max;
  vector<Cell*> remainingDaughters(mTree[0].begin()->second.daughterPtrs), newDaughters;
  vector<int> result;

  // Determine the grid cell range on the master level.
  const CellKey ix_master_min = ix_master == 0        ? ix_master : ix_master - 1;
  const CellKey iy_master_min = iy_master == 0        ? iy_master : iy_master - 1;
  const CellKey iz_master_min = iz_master == 0        ? iz_master : iz_master - 1;
  const CellKey ix_master_max = ix_master == max1dKey ? ix_master : ix_master + 1;
  const CellKey iy_master_max = iy_master == max1dKey ? iy_master : iy_master + 1;
  const CellKey iz_master_max = iz_master == max1dKey ? iz_master : iz_master + 1;

  // Walk the tree until we run out of daughters to check.
  while (remainingDaughters.size() > 0) {
    newDaughters = vector<Cell*>();
    ++ilevel;

    // Find the target range of keys on this level.
    ix_min = shiftKeyLevel(ix_master_min, masterLevel, ilevel);
    iy_min = shiftKeyLevel(iy_master_min, masterLevel, ilevel);
    iz_min = shiftKeyLevel(iz_master_min, masterLevel, ilevel);
    ix_max = shiftKeyLevel(ix_master_max, masterLevel, ilevel);
    iy_max = shiftKeyLevel(iy_master_max, masterLevel, ilevel);
    iz_max = shiftKeyLevel(iz_master_max, masterLevel, ilevel);
    CHECK(ix_min <= ix_max and ix_max <= max1dKey);
    CHECK(iy_min <= iy_max and iy_max <= max1dKey);
    CHECK(iz_min <= iz_max and iz_max <= max1dKey);
    
    // Walk the candidate daughters on this level.
    for (typename vector<Cell*>::const_iterator itr = remainingDaughters.begin();
         itr != remainingDaughters.end();
         ++itr) {
      const Cell& cell = **itr;
      
      // Is this daughter in range?
      if (keyInRange(cell.key, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max)) {
        
        // Copy this cells members to the result.
        copy(cell.members.begin(), cell.members.end(), back_inserter(result));

        // Add any daughters of this cell to our candidates to check on the next level.
        copy(cell.daughterPtrs.begin(), cell.daughterPtrs.end(), back_inserter(newDaughters));
      }
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Shift a 1D key from one level to another.
//------------------------------------------------------------------------------
template<typename Dimension>
typename TreeNeighbor<Dimension>::CellKey
TreeNeighbor<Dimension>::
shiftKeyLevel(const typename TreeNeighbor<Dimension>::CellKey& ix,
              const typename TreeNeighbor<Dimension>::LevelKey& level0,
              const typename TreeNeighbor<Dimension>::LevelKey& level1) const {
  return (level1 <= level0 ? 
          ix >> (level0 - level1) :
          ix << (level1 - level0));
}

//------------------------------------------------------------------------------
// Check if the given key is in a specified range.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TreeNeighbor<Dimension>::
keyInRange(const typename TreeNeighbor<Dimension>::CellKey& key,
           const typename TreeNeighbor<Dimension>::CellKey& ix_min,
           const typename TreeNeighbor<Dimension>::CellKey& iy_min,
           const typename TreeNeighbor<Dimension>::CellKey& iz_min,
           const typename TreeNeighbor<Dimension>::CellKey& ix_max,
           const typename TreeNeighbor<Dimension>::CellKey& iy_max,
           const typename TreeNeighbor<Dimension>::CellKey& iz_max) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  return (ix >= ix_min and ix <= ix_max and
          iy >= iy_min and iy <= iy_max and
          iz >= iz_min and iz <= iz_max);
}

//------------------------------------------------------------------------------
// Define our static members.
//------------------------------------------------------------------------------
template<typename Dimension> const unsigned TreeNeighbor<Dimension>::num1dbits = 21U;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::max1dKey = 1U << TreeNeighbor<Dimension>::num1dbits;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::xkeymask = (1U << TreeNeighbor<Dimension>::num1dbits) - 1U;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::ykeymask = TreeNeighbor<Dimension>::xkeymask << TreeNeighbor<Dimension>::num1dbits;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::zkeymask = TreeNeighbor<Dimension>::ykeymask << TreeNeighbor<Dimension>::num1dbits;

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TreeNeighbor<Dim<1> >;

}
}
