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
#include "Boundary/mapPositionThroughPlanes.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace NeighborSpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Compute the vertex coordinates for a cell.
//------------------------------------------------------------------------------
// 1D
inline
vector<Dim<1>::Vector>
findCellVertices(const Dim<1>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& iy,
                 const uint64_t& iz) {
  typedef Dim<1>::Vector Vector;
  const double cellSize = boxLength/(1U << ilevel);
  vector<Vector> result;
  result.push_back(xmin + Vector(ix       *cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize));
  return result;
}

// 2D
inline
vector<Dim<2>::Vector>
findCellVertices(const Dim<2>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& iy,
                 const uint64_t& iz) {
  typedef Dim<2>::Vector Vector;
  const double cellSize = boxLength/(1U << ilevel);
  vector<Vector> result;
  result.push_back(xmin + Vector(ix       *cellSize, iy       *cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, iy       *cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize));
  result.push_back(xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize));
  return result;
}

// 3D
inline
vector<Dim<3>::Vector>
findCellVertices(const Dim<3>::Vector& xmin,
                 const double& boxLength,
                 const uint32_t& ilevel,
                 const uint64_t& ix,
                 const uint64_t& iy,
                 const uint64_t& iz) {
  typedef Dim<3>::Vector Vector;
  const double cellSize = boxLength/(1U << ilevel);
  vector<Vector> result;
  result.push_back(xmin + Vector(ix       *cellSize, iy       *cellSize, iz       *cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, iz       *cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, iz       *cellSize));
  result.push_back(xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, iz       *cellSize));
  result.push_back(xmin + Vector(ix       *cellSize, iy       *cellSize, (iz + 1U)*cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, (iz + 1U)*cellSize));
  result.push_back(xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize));
  result.push_back(xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize));
  return result;
}

//------------------------------------------------------------------------------
// Squeeze the vertices of a cell a bit closer together.
//------------------------------------------------------------------------------
// 1D
inline
void
squeezeCell(vector<Dim<1>::Vector>& vertices,
            const double& smidgen) {
  typedef Dim<1>::Vector Vector;
  REQUIRE(vertices.size() == 2);
  const Vector centroid = 0.5*smidgen*(vertices[0] + vertices[1]);
  vertices[0] = centroid + (1.0 - smidgen)*vertices[0];
  vertices[1] = centroid + (1.0 - smidgen)*vertices[1];
}

// 2D
inline
void
squeezeCell(vector<Dim<2>::Vector>& vertices,
            const double& smidgen) {
  typedef Dim<2>::Vector Vector;
  REQUIRE(vertices.size() == 4);
  const Vector centroid = 0.25*smidgen*(vertices[0] + vertices[1] + 
                                        vertices[2] + vertices[3]);
  vertices[0] = centroid + (1.0 - smidgen)*vertices[0];
  vertices[1] = centroid + (1.0 - smidgen)*vertices[1];
  vertices[2] = centroid + (1.0 - smidgen)*vertices[2];
  vertices[3] = centroid + (1.0 - smidgen)*vertices[3];
}

// 3D
inline
void
squeezeCell(vector<Dim<3>::Vector>& vertices,
            const double& smidgen) {
  typedef Dim<3>::Vector Vector;
  REQUIRE(vertices.size() == 8);
  const Vector centroid = 0.125*smidgen*(vertices[0] + vertices[1] + vertices[2] + vertices[3] +
                                         vertices[4] + vertices[5] + vertices[6] + vertices[7]);
  vertices[0] = centroid + (1.0 - smidgen)*vertices[0];
  vertices[1] = centroid + (1.0 - smidgen)*vertices[1];
  vertices[2] = centroid + (1.0 - smidgen)*vertices[2];
  vertices[3] = centroid + (1.0 - smidgen)*vertices[3];
  vertices[4] = centroid + (1.0 - smidgen)*vertices[4];
  vertices[5] = centroid + (1.0 - smidgen)*vertices[5];
  vertices[6] = centroid + (1.0 - smidgen)*vertices[6];
  vertices[7] = centroid + (1.0 - smidgen)*vertices[7];
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeNeighbor<Dimension>::
TreeNeighbor(NodeList<Dimension>& nodeList,
             const NeighborSearchType searchType,
             const double kernelExtent,
             const Vector& xmin,
             const Vector& xmax):
  Neighbor<Dimension>(nodeList, searchType, kernelExtent),
  mBoxLength((xmax - xmin).maxElement()),
  mGridLevelConst0(log(mBoxLength/kernelExtent)/log(2.0)),
  mXmin(xmin),
  mXmax(xmax),
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

  // Get the master and coarse lists.
  vector<int>& masterList = this->accessMasterList();
  vector<int>& coarseList = this->accessCoarseNeighborList();
  masterList = vector<int>();
  coarseList = vector<int>();
  if (mTree.size() > 0) {
    CHECK(mTree[0].size() == 1);
    CHECK(mTree[0].begin()->second.members.size() == 0);

    // Declare a bunch of variables we're going to need.
    LevelKey ilevel = 0;
    CellKey ix, iy, iz;
    double cellSize;
    vector<Cell*> remainingDaughters(mTree[0].begin()->second.daughterPtrs), newDaughters;

    // Walk the tree, looking for any master cells that are in range of the 
    // entrance plane.
    while (remainingDaughters.size() > 0) {
      newDaughters = vector<Cell*>();
      ++ilevel;
      cellSize = mBoxLength/(1U << ilevel);
    
      // Walk the candidates.
      for (typename vector<Cell*>::const_iterator itr = remainingDaughters.begin();
           itr != remainingDaughters.end();
           ++itr) {
        const Cell& cell = **itr;

        // Check if we're in range of either plane.
        const bool entranceCheck = (this->distanceToCell(ilevel, cell.key, enterPlane) <= cellSize);
        const bool exitCheck = (this->distanceToCell(ilevel, cell.key, exitPlane) <= cellSize);

        // If so, add the daughters to check on the next pass.
        if (entranceCheck or exitCheck) copy(cell.daughterPtrs.begin(), cell.daughterPtrs.end(), back_inserter(newDaughters));

        // Does this cell have members in range of the entrance plane?
        if (entranceCheck and cell.members.size() > 0) {
          copy(cell.members.begin(), cell.members.end(), back_inserter(masterList));

          // Map the cell key through to the exit plane (which may result in more
          // than one equivalent key).  Then find all the neighbors for those cells,
          // and add them to the coarse set.
          const vector<CellKey> mappedKeys = this->mapKey(ilevel, cell.key, enterPlane, exitPlane);
          for (typename vector<CellKey>::const_iterator keyItr = mappedKeys.begin();
               keyItr != mappedKeys.end();
               ++keyItr) {
            this->extractCellIndices(*keyItr, ix, iy, iz);
            const vector<int> neighbors = this->findTreeNeighbors(ilevel, ix, iy, iz);
            copy(neighbors.begin(), neighbors.end(), back_inserter(coarseList));
          }
        }

        // Does this cell have members in range of the exit plane?
        if (exitCheck and cell.members.size() > 0) {
          copy(cell.members.begin(), cell.members.end(), back_inserter(coarseList));

          // Map the cell key through to the entrance plane.  Any nodes we interact
          // with on that side are potential masters.
          const vector<CellKey> mappedKeys = this->mapKey(ilevel, cell.key, exitPlane, enterPlane);
          for (typename vector<CellKey>::const_iterator keyItr = mappedKeys.begin();
               keyItr != mappedKeys.end();
               ++keyItr) {
            this->extractCellIndices(*keyItr, ix, iy, iz);
            const vector<int> masters = this->findTreeNeighbors(ilevel, ix, iy, iz);
            copy(masters.begin(), masters.end(), back_inserter(masterList));
            // TreeLevel::const_iterator cellItr = mTree[ilevel].find(*keyItr);
            // if (cellItr != mTree[ilevel].end()) {
            //   this->extractCellIndices(*keyItr, ix, iy, iz);
            //   const vector<int> masters = this->findTreeNeighbors(ilevel, ix, iy, iz);
            //   copy(masters.begin(), masters.end(), back_inserter(masterList));
            // }
          }
        }
      }

      // Update the set of daughters to check on the next pass.
      remainingDaughters = newDaughters;
    }

    // Remove duplicates from the master and coarse sets.
    sort(masterList.begin(), masterList.end());
    masterList.erase(unique(masterList.begin(), masterList.end()), masterList.end());
    sort(coarseList.begin(), coarseList.end());
    coarseList.erase(unique(coarseList.begin(), coarseList.end()), coarseList.end());

    // Ghost have to be allowed to be master for boundary conditions to work!
    // // We don't allow ghost nodes to be masters.
    // const int firstGhostNode = this->nodeList().firstGhostNode();
    // masterList.erase(lower_bound(masterList.begin(), masterList.end(), firstGhostNode), masterList.end());
  }
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

  // // Recompute the current box size.  We assume xmin & xmax have
  // // already been set.
  // mBoxLength = (this->mXmax - this->mXmin).maxElement();
  // CHECK(mBoxLength > 0.0);
  // mGridLevelConst0 = log(mBoxLength/this->kernelExtent())/log(2.0);

  // Walk all the nodes and add them to the tree.
  const size_t n = nodes.numNodes();
  for (unsigned i = 0; i != n; ++i) {
    this->addNodeToTree(positions(i), H(i), i);
  }

  // Set the daughter pointers.
  constructDaughterPtrs(mTree);

  // Force the node extents to be calculated.
  this->setNodeExtents();

  ENSURE(valid());
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
  REQUIRE2(this->kernelExtent()*h <= mBoxLength,
           "h larger than box size: " << this->kernelExtent()*h << " " << mBoxLength);
  const unsigned result = std::max(0, 
                                   std::min(int(num1dbits) - 1,
                                            int(mGridLevelConst0 - log(h)/log(2.0))));
  ENSURE(result < num1dbits);
  ENSURE2((result == 0) or (result == num1dbits - 1) or
          (fuzzyLessThanOrEqual(h*this->kernelExtent(), mBoxLength/(1U << result), 1.0e-10) and
           fuzzyGreaterThanOrEqual(h*this->kernelExtent(), mBoxLength/(1U << (result + 1U)), 1.0e-10)),
          result << " " << h*this->kernelExtent() << " in? ["
          << mBoxLength/(1U << (result + 1U)) << " " 
          << mBoxLength/(1U << result) << "]\n"
          << mBoxLength << " " << mGridLevelConst0);
          
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
  if (globalTree) nlevels = allReduce(nlevels, MPI_MAX, Communicator::communicator());

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
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
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
  if (globalTree) nlevels = allReduce(nlevels, MPI_MAX, Communicator::communicator());

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
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
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
// xmin
//------------------------------------------------------------------------------
template<typename Dimension>
const typename Dimension::Vector&
TreeNeighbor<Dimension>::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
template<typename Dimension>
const typename Dimension::Vector&
TreeNeighbor<Dimension>::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// boxlength
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeNeighbor<Dimension>::
boxLength() const {
  return mBoxLength;
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
  packElement(mXmin, buffer);
  packElement(mXmax, buffer);
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
  unpackElement(mXmin, bufItr, endItr);
  unpackElement(mXmax, bufItr, endItr);
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
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
buildCellKey(const typename TreeNeighbor<Dimension>::LevelKey ilevel,
             const typename TreeNeighbor<Dimension>::Vector& xi,
             typename TreeNeighbor<Dimension>::CellKey& key,
             typename TreeNeighbor<Dimension>::CellKey& ix,
             typename TreeNeighbor<Dimension>::CellKey& iy,
             typename TreeNeighbor<Dimension>::CellKey& iz) const {
//   REQUIRE2(xi.x() >= this->mXmin.x() and xi.x() <= this->mXmax.x(), xi << " " << this->mXmin << " " << this->mXmax);
//   REQUIRE2(xi.y() >= this->mXmin.y() and xi.y() <= this->mXmax.y(), xi << " " << this->mXmin << " " << this->mXmax);
//   REQUIRE2(xi.z() >= this->mXmin.z() and xi.z() <= this->mXmax.z(), xi << " " << this->mXmin << " " << this->mXmax);
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (xi.x() - this->mXmin.x())/mBoxLength)) * ncell));
  iy = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (xi.y() - this->mXmin.y())/mBoxLength)) * ncell));
  iz = std::min(maxcell, CellKey(std::max(0.0, std::min(1.0, (xi.z() - this->mXmin.z())/mBoxLength)) * ncell));
  key = ((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
extractCellIndices(const typename TreeNeighbor<Dimension>::CellKey& key,
                   typename TreeNeighbor<Dimension>::CellKey& ix,
                   typename TreeNeighbor<Dimension>::CellKey& iy,
                   typename TreeNeighbor<Dimension>::CellKey& iz) const {
  ix = key & xkeymask;
  iy = (key & ykeymask) >> num1dbits;
  iz = (key & zkeymask) >> 2*num1dbits;
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
addDaughter(typename TreeNeighbor<Dimension>::Cell& cell,
            const typename TreeNeighbor<Dimension>::CellKey daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= (1U << Dimension::nDim));
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
addNodeToTree(const typename Dimension::Vector& xi,
              const typename Dimension::SymTensor& Hi,
              const unsigned i) {
  mTree.reserve(num1dbits); // This is necessary to avoid memory errors!

  // Determine the level for this point.
  const LevelKey homeLevel = this->gridLevel(Hi);

  LevelKey ilevel = 0;
  CellKey key, parentKey, otherKey, ix, iy, iz;
  typename TreeLevel::iterator itr;

  // First walk the tree for all levels above our native level.
  while (ilevel <= homeLevel) {

    // Do we need to add another level to the tree?
    if (ilevel == mTree.size()) mTree.push_back(TreeLevel());

    // Create the key for the cell containing this particle on this level.
    buildCellKey(ilevel, xi, key, ix, iy, iz);
    itr = mTree[ilevel].find(key);

    // Is this a new cell?
    if (itr == mTree[ilevel].end()) {
      mTree[ilevel][key] = Cell(key);
      itr = mTree[ilevel].find(key);
    }

    // Link this cell as a daughter of its parent.
    if (ilevel > 0) {
      CHECK(mTree[ilevel - 1].find(parentKey) != mTree[ilevel - 1].end());
      addDaughter(mTree[ilevel - 1][parentKey], key);
    }

    // Is this the final level for this node?
    if (ilevel == homeLevel) {
      itr->second.members.push_back(i);
    }

    // Prepare for the next level.
    parentKey = key;
    ++ilevel;
  }

}

//------------------------------------------------------------------------------
// Construct the daughter pointers in a tree.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
constructDaughterPtrs(typename TreeNeighbor<Dimension>::Tree& tree) const {
  const unsigned nlevels = tree.size();
  const unsigned n = nlevels > 0 ? nlevels - 1 : nlevels;
  unsigned ilevel, ilevel1;
  for (ilevel = 0; ilevel != n; ++ilevel) {
    ilevel1 = ilevel + 1;
    for (typename TreeLevel::iterator itr = tree[ilevel].begin();
         itr != tree[ilevel].end();
         ++itr) {
      Cell& cell = itr->second;
      cell.daughterPtrs = std::vector<Cell*>();
      for (typename std::vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) {
        cell.daughterPtrs.push_back(&(tree[ilevel1][*ditr]));
      }
      CHECK(cell.daughters.size() == cell.daughterPtrs.size());
    }
  }
}

//------------------------------------------------------------------------------
// Set the master & coarse neighbor sets by walking the tree.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeNeighbor<Dimension>::
setTreeMasterList(const typename Dimension::Vector& position,
                  const double& h) {

  // Set the working master grid level and cell.
  CellKey masterKey, ix_master, iy_master, iz_master;
  const LevelKey masterLevel = this->gridLevel(h);
  buildCellKey(masterLevel, position, masterKey, ix_master, iy_master, iz_master);
  CHECK(masterLevel >= 0 and masterLevel < num1dbits);

  // Grab the lists we're going to fill in.
  vector<int>& masterList = this->accessMasterList();
  vector<int>& coarseNeighborList = this->accessCoarseNeighborList();
  masterList = vector<int>();
  coarseNeighborList = vector<int>();

  // Set the master list.
  if (mTree.size() > masterLevel) {
    typename TreeLevel::const_iterator masterItr = mTree[masterLevel].find(masterKey);
    masterList = (masterItr == mTree[masterLevel].end() ? vector<int>() : masterItr->second.members);
  }

  // Set the coarse list.
  if (mTree.size() > 0) {
    coarseNeighborList = this->findTreeNeighbors(masterLevel, ix_master, iy_master, iz_master);
  }

  // Remove all ghost nodes from the master list.
  const unsigned firstGhostNode = this->nodeList().firstGhostNode();
  sort(masterList.begin(), masterList.end());
  masterList.erase(lower_bound(masterList.begin(), masterList.end(), firstGhostNode), masterList.end());

  // Post conditions.
  ENSURE(coarseNeighborList.size() >= this->masterList().size());
  ENSURE(masterList.size() == 0 or *max_element(masterList.begin(), masterList.end()) < firstGhostNode);
}

//------------------------------------------------------------------------------
// Produce the refined list of potential neighbors for a single node.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeNeighbor<Dimension>::
setTreeRefineNeighborList(const typename Dimension::Vector& position,
                          const typename Dimension::SymTensor& H) {

  // // Determine the maximum extent of this H tensor in each dimension.
  // const Vector extent = this->HExtent(H, this->kernelExtent());
  // const Vector minExtent = position - extent;
  // const Vector maxExtent = position + extent;

  // Use precull to set the refined neighbor list.
  const std::vector<int>& coarseList = this->coarseNeighborList();
  std::vector<int>& refineList = this->accessRefineNeighborList();
  refineList = coarseList; // this->precullList(position, position, minExtent, maxExtent, coarseList);
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

  // Declare variables.
  LevelKey ilevel = 0;
  CellKey ix, iy, iz, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max, delta;
  vector<Cell*> remainingDaughters(mTree[0].begin()->second.daughterPtrs), newDaughters;
  vector<int> result;

  // Walk the tree until we run out of daughters to check.
  CHECK2(mTree[0].begin()->second.members.size() == 0, mTree[0].begin()->second.members.size());
  while (remainingDaughters.size() > 0) {
    newDaughters = vector<Cell*>();
    ++ilevel;
    delta = (ilevel <= masterLevel ? 1U : (1U << (ilevel - masterLevel)));

    // Find the target range of keys on this level.
    ix = shiftKeyLevel(ix_master, masterLevel, ilevel);
    iy = shiftKeyLevel(iy_master, masterLevel, ilevel);
    iz = shiftKeyLevel(iz_master, masterLevel, ilevel);
    ix_min = (ix > delta              ? ix - delta : 0U);
    iy_min = (iy > delta              ? iy - delta : 0U);
    iz_min = (iz > delta              ? iz - delta : 0U);
    ix_max = ((max1dKey - ix) > delta ? ix + 2*delta - 1U : max1dKey);
    iy_max = ((max1dKey - iy) > delta ? iy + 2*delta - 1U : max1dKey);
    iz_max = ((max1dKey - iz) > delta ? iz + 2*delta - 1U : max1dKey);
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

    // Update the daughters to check on the next pass.
    remainingDaughters = newDaughters;
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
// Find the minimum distance from a cell to a plane.
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeNeighbor<Dimension>::
distanceToCell(const typename TreeNeighbor<Dimension>::LevelKey& ilevel,
               const typename TreeNeighbor<Dimension>::CellKey& key,
               const GeomPlane<Dimension>& plane) const {
  CellKey ix, iy, iz;
  this->extractCellIndices(key, ix, iy, iz);
  const vector<Vector> cellVertices = findCellVertices(mXmin, mBoxLength, 
                                                       ilevel, ix, iy, iz);
  const unsigned n = cellVertices.size();
  CHECK(n == (1U << Dimension::nDim));
  unsigned i;
  double dist, minDist = plane.signedDistance(cellVertices[0]);
  for (unsigned i = 1; i < cellVertices.size(); ++i) {
    dist = plane.signedDistance(cellVertices[i]);
    if (dist*minDist < 0.0) {
      minDist = 0.0;
    } else if (abs(dist) < abs(minDist)) {
      minDist = dist;
    }
  }
  return abs(minDist);
}

//------------------------------------------------------------------------------
// Map a cell key/level through an entrance/exit plane pair, returning the set 
// of cells that overlap the mapped cell.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<typename TreeNeighbor<Dimension>::CellKey>
TreeNeighbor<Dimension>::
mapKey(const typename TreeNeighbor<Dimension>::LevelKey& ilevel,
       const typename TreeNeighbor<Dimension>::CellKey& key,
       const GeomPlane<Dimension>& enterPlane,
       const GeomPlane<Dimension>& exitPlane) const {
  CellKey ix, iy, iz;
  this->extractCellIndices(key, ix, iy, iz);
  vector<Vector> vertices = findCellVertices(mXmin, mBoxLength,
                                             ilevel, ix, iy, iz);
  const unsigned n = vertices.size();
  CHECK(n == (1U << Dimension::nDim));

  // Find the range of (ix,iy,iz) the mapped vertices cover.
  CellKey ixmin = max1dKey, 
          iymin = max1dKey, 
          izmin = max1dKey,
          ixmax = CellKey(0),
          iymax = CellKey(0),
          izmax = CellKey(0),
          newKey;
  unsigned i;
  for (i = 0; i != n; ++i) {
    buildCellKey(ilevel,
                 mapPositionThroughPlanes(vertices[i], enterPlane, exitPlane),
                 newKey, ix, iy, iz);
    ixmin = min(ixmin, ix);
    iymin = min(iymin, iy);
    izmin = min(izmin, iz);
    ixmax = max(ixmax, ix);
    iymax = max(iymax, iy);
    izmax = max(izmax, iz);
  }

  // Now fill in the set of cells we map to.
  vector<CellKey> result;
  for (ix = ixmin; ix <= ixmax; ++ix) {
    for (iy = iymin; iy <= iymax; ++iy) {
      for (iz = izmin; iz <= izmax; ++iz) {
        result.push_back((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
                         (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
                         (std::max(CellKey(0), std::min(max1dKey, ix))));
      }
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Internal consistency checks.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TreeNeighbor<Dimension>::
valid() const {
  bool result = true;

  // Make sure each node is listed only once.
  map<unsigned, unsigned> nodeCount;
  for (typename Tree::const_iterator levelItr = mTree.begin(); 
       levelItr != mTree.end();
       ++levelItr) {
    for (typename TreeLevel::const_iterator cellItr = levelItr->begin();
         cellItr != levelItr->end();
         ++cellItr) {
      const Cell& cell = cellItr->second;
      for (vector<int>::const_iterator iitr = cell.members.begin();
           iitr != cell.members.end();
           ++iitr) {
        map<unsigned, unsigned>::iterator itr = nodeCount.find(*iitr);
        if (itr == nodeCount.end()) {
          itr->second = 1;
        } else {
          ++(itr->second);
        }
      }
    }
  }

  for (typename map<unsigned, unsigned>::const_iterator itr = nodeCount.begin();
       itr != nodeCount.end();
       ++itr) {
    const unsigned i = itr->first;
    const unsigned count = itr->second;
    if (count != 1) {
      cerr << "TreeNeighor::valid failing test of nodes uniquely assigned to cell: " << i << " " << count << endl;
      result = false;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Define our static members.
//------------------------------------------------------------------------------
template<typename Dimension> const unsigned TreeNeighbor<Dimension>::num1dbits = 21U;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::max1dKey = 1U << TreeNeighbor<Dimension>::num1dbits;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::xkeymask = (1U << TreeNeighbor<Dimension>::num1dbits) - 1U;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::ykeymask = TreeNeighbor<Dimension>::xkeymask << TreeNeighbor<Dimension>::num1dbits;
template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::zkeymask = TreeNeighbor<Dimension>::ykeymask << TreeNeighbor<Dimension>::num1dbits;

}
}
