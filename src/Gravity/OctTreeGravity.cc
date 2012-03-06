//---------------------------------Spheral++----------------------------------//
// OctTreeGravity -- An implementation of the tree n-body gravity solver.
//
// Created by JMO, 2012-02-28
//----------------------------------------------------------------------------//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "OctTreeGravity.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/PairComparisons.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "DBC.hh"

namespace Spheral {
namespace GravitySpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
OctTreeGravity::
OctTreeGravity(const double G,
               const double softeningLength,
               const double opening,
               const double ftimestep):
  mG(G),
  mSofteningLength(softeningLength),
  mOpening(opening),
  mftimestep(ftimestep),
  mBoxLength(0.0),
  mMaxCellDensity(0.0),
  mXmin(),
  mXmax(),
  mTree(),
  mPotential(FieldList<Dim<3>, Scalar>::Copy),
  mExtraEnergy(0.0) {
  VERIFY(G > 0.0);
  VERIFY(opening > 0.0);
  VERIFY(softeningLength > 0.0);
  VERIFY(ftimestep > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
OctTreeGravity::
~OctTreeGravity() {
}

//------------------------------------------------------------------------------
// Register some extra state.
//------------------------------------------------------------------------------
void
OctTreeGravity::
registerState(DataBase<Dim<3> >& dataBase,
              State<Dim<3> >& state) {
  GenericBodyForce<Dim<3> >::registerState(dataBase, state);
  state.enrollFieldList(mPotential);
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
void 
OctTreeGravity::
evaluateDerivatives(const Dim<3>::Scalar time,
                    const Dim<3>::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dim<3> >& state,
                    StateDerivatives<Dim<3> >& derivs) const {

  // Access the pertinent fields in the database.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);

  // Get the acceleration and position change vectors we'll be modifying.
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;
  mPotential = 0.0;

#ifdef USE_MPI

  // Get the processor information.
  const unsigned rank = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();

  // Pack up the local tree.
  vector<char> localBuffer, buffer;
  this->serialize(mTree, localBuffer);

  // Walk each process, and have them send their local tree info to all other
  // domains.
  unsigned bufSize;
  vector<char>::const_iterator bufItr;
  Tree tree;
  for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {

    // Broadcast the send processor's tree.
    buffer = localBuffer;
    bufSize = buffer.size();
    MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
    buffer.resize(bufSize);
    MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
    tree = Tree();
    bufItr = buffer.begin();
    this->deserialize(tree, bufItr, buffer.end());
    CHECK(bufItr == buffer.end());

    // Add this partial tree's contribution to our local forces and potentials.
    applyTreeForces(tree, mass, position, DxDt, DvDt, mPotential);
  }

#else

  // Apply the forces from our local tree.
  applyTreeForces(mTree, mass, position, DxDt, DvDt, mPotential);

#endif

  // Set the motion to be Lagrangian.
  DxDt.assignFields(velocity);
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
void 
OctTreeGravity::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}

//------------------------------------------------------------------------------
// Intialize the package before evaluateDerivatives is called.
// For OctTreeGravity, this is where we build the current tree.
//------------------------------------------------------------------------------
void 
OctTreeGravity::
initialize(const Scalar time,
           const Scalar dt,
           const DataBase<Dimension>& db,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // For now we're not going to be clever about trying to patch an existing tree,
  // but instead we'll build it from scratch every time.
  mTree = Tree();

  // Access to pertinent fields in the database.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const size_t numNodeLists = mass.numFields();

  // Determine the box size.
  globalBoundingBox(position, mXmin, mXmax, false);
  CHECK(mXmin.x() <= mXmax.x());
  CHECK(mXmin.y() <= mXmax.y());
  CHECK(mXmin.z() <= mXmax.z());
  mBoxLength = (mXmax - mXmin).maxAbsElement();
  CHECK(mBoxLength > 0.0);

  // Walk all the internal nodes and add them to a temporary MapTree.
  MapTree tree;
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const size_t n = mass[nodeListi]->numInternalElements();
    for (size_t i = 0; i != n; ++i) {
      this->addNodeToTree(tree, mass(nodeListi, i), position(nodeListi, i));
    }
  }

  // Copy the result to the internal tree, and sort it by the TreeKeys.
  for (MapTree::const_iterator itr = tree.begin();
       itr != tree.end();
       ++itr) mTree.push_back(*itr);
  sort(mTree.begin(), mTree.end(), TreeComparatorLessThan());

#ifdef USE_MPI
  {

    // We require the true global center of mass for each cell, which means
    // in parallel we need to exchange the local trees and build up this information.
    // We also build the global total mass in each cell at the same time.
    // Get the processor information.
    const unsigned rank = Process::getRank();
    const unsigned numProcs = Process::getTotalNumberOfProcesses();

    // Pack up the local tree.
    vector<char> localBuffer, buffer;
    this->serialize(mTree, localBuffer);

    // Walk each process, and have them send their local tree info to all other
    // domains.
    unsigned bufSize;
    vector<char>::const_iterator bufItr;
    Tree tree;
    for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {

      // Broadcast the send processor's tree.
      buffer = localBuffer;
      bufSize = buffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
      tree = Tree();
      bufItr = buffer.begin();
      this->deserialize(tree, bufItr, buffer.end());
      CHECK(bufItr == buffer.end());

      // Augment any of our local cell calculations for mass and center of mass.
      for (Tree::iterator myItr = mTree.begin();
           myItr != mTree.end();
           ++myItr) {
        const TreeKey key = myItr->first;
        const Tree::const_iterator otherItr = lower_bound(tree.begin(),
                                                          tree.end(),
                                                          key,
                                                          TreeComparatorLessThan());
        if (otherItr != tree.end() and otherItr->first == key) {
          Cell& cell = myItr->second;
          const Cell& otherCell = otherItr->second;
          cell.xcm = (cell.Mglobal*cell.xcm + otherCell.M*otherCell.xcm)/(cell.Mglobal + otherCell.M);
          cell.Mglobal += otherCell.M;
        }
      }
    }

  }
#endif

  // Make a final pass over the cells and fill in the distance between the 
  // center of mass and the geometric center.
  // We also squirrel away the maximum effective cell density for our timestep
  // determination.
  unsigned ilevel;
  CellKey ckey, ix, iy, iz;
  double cellsize, cellvol;
  mMaxCellDensity = 0.0;
  for (Tree::iterator itr = mTree.begin();
       itr != mTree.end();
       ++itr) {
    ilevel = itr->first.first;
    ckey = itr->first.second;
    Cell& cell = itr->second;
    cellsize = mBoxLength/(1U << ilevel);
    cellvol = cellsize*cellsize*cellsize;
    extractCellIndices(ckey, ix, iy, iz);

    // Update the distance between the cell's center of mass and geometric center.
    cell.rcm2cc = (cell.xcm - (mXmin + Vector((ix + 0.5)*cellsize,
                                              (iy + 0.5)*cellsize,
                                              (iz + 0.5)*cellsize))).magnitude();
    CHECK(cell.rcm2cc < 1.74*cellsize);

    // Update the maximum effective cell density.
    mMaxCellDensity = max(mMaxCellDensity, cell.Mglobal/cellvol);
  }

  // Get the global max cell density.
  mMaxCellDensity = allReduce(mMaxCellDensity, MPI_MAX, MPI_COMM_WORLD);
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
OctTreeGravity::TimeStepType
OctTreeGravity::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  REQUIRE(mMaxCellDensity > 0.0);

  // We use the gravitational dynamical time (sqrt(G/rho)) to estimate the 
  // necessary timestep.
  const double dt = mftimestep * sqrt(1.0/(mG*mMaxCellDensity));

  stringstream reasonStream;
  reasonStream << "OctTreeGravity: sqrt(1/(G rho)) = sqrt(1/("
               << mG << " * " << mMaxCellDensity
               << ")) = " << dt << ends;
  return TimeStepType(dt, reasonStream.str());
}

//------------------------------------------------------------------------------
// extraEnergy
//------------------------------------------------------------------------------
OctTreeGravity::Scalar 
OctTreeGravity::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
const FieldList<Dim<3>, OctTreeGravity::Scalar>&
OctTreeGravity::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
std::string
OctTreeGravity::
dumpTree() const {

  stringstream ss;
  CellKey key, ix, iy, iz;
  Tree::const_iterator itr = mTree.begin(), endItr = mTree.end(), lastItr;
  const unsigned nlevels = mTree.size() > 0 ? mTree.back().first.first + 1 : 0;

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    lastItr = lower_bound(itr, endItr, TreeKey(ilevel+1, 0U), TreeComparatorLessThan());
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " 
       << distance(itr, lastItr) << "\n";
    for (; itr != lastItr; ++itr) {
      key = itr->first.second;
      const Cell& cell = itr->second;
      extractCellIndices(key, ix, iy, iz);
      ss << "    Cell key=" << key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << ")\n"
         << "         xcm=" << cell.xcm << " rcm2cc=" << cell.rcm2cc << " M=" << cell.M << "\n"
         << "         daughters = ( ";
      for (vector<TreeKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) ss << ditr->second << " ";
      ss << ")\n"
         << "         nodes = [";
      for (unsigned k = 0; k != cell.masses.size(); ++k) ss << " ("
                                                            << cell.masses[k] << " "
                                                            << cell.positions[k] << ")";
      ss <<" ]\n";
    }
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
std::string
OctTreeGravity::
dumpTreeStatistics() const {

  stringstream ss;
  CellKey key, ix, iy, iz;
  Tree::const_iterator itr = mTree.begin(), endItr = mTree.end(), lastItr;
  const unsigned nlevels = mTree.size() > 0 ? mTree.back().first.first + 1 : 0;
  unsigned ncells, nparticles;

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    lastItr = lower_bound(itr, endItr, TreeKey(ilevel+1, 0U), TreeComparatorLessThan());
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " 
       << distance(itr, lastItr) << "\n";
    ncells = 0;
    nparticles = 0;
    for (; itr != lastItr; ++itr) {
      key = itr->first.second;
      const Cell& cell = itr->second;
      ++ncells;
      nparticles += cell.masses.size();
    }
    ss << "         : ncells = " << ncells << "\n"
       << "         : nparts = " << nparticles << "\n";
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
double
OctTreeGravity::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// opening
//------------------------------------------------------------------------------
double
OctTreeGravity::
opening() const {
  return mOpening;
}

void
OctTreeGravity::
opening(const double x) {
  VERIFY(x > 0.0);
  mOpening = x;
}

//------------------------------------------------------------------------------
// softeningLength
//------------------------------------------------------------------------------
double
OctTreeGravity::
softeningLength() const {
  return mSofteningLength;
}

void
OctTreeGravity::
softeningLength(const double x) {
  VERIFY(x > 0.0);
  mSofteningLength = x;
}

//------------------------------------------------------------------------------
// ftimestep
//------------------------------------------------------------------------------
double
OctTreeGravity::
ftimestep() const {
  return mftimestep;
}

void
OctTreeGravity::
ftimestep(const double x) {
  VERIFY(x > 0.0);
  mftimestep = x;
}

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
OctTreeGravity::Vector
OctTreeGravity::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
OctTreeGravity::Vector
OctTreeGravity::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// maxCellDensity
//------------------------------------------------------------------------------
double
OctTreeGravity::
maxCellDensity() const {
  return mMaxCellDensity;
}

//------------------------------------------------------------------------------
// Apply the forces from a tree.
//------------------------------------------------------------------------------
void 
OctTreeGravity::
applyTreeForces(const Tree& tree,
                const FieldSpace::FieldList<Dimension, Scalar>& mass,
                const FieldSpace::FieldList<Dimension, Vector>& position,
                FieldSpace::FieldList<Dimension, Vector>& DxDt,
                FieldSpace::FieldList<Dimension, Vector>& DvDt,
                FieldSpace::FieldList<Dimension, Scalar>& potential) const {

  const double soft2 = mSofteningLength*mSofteningLength;
  const unsigned numNodeLists = mass.numFields();
  const unsigned numNodes = mass.numInternalNodes();

  // We'll always be starting with the daughters of the root level.
  const Cell& rootCell = tree[0].second;

  // Declare stuff we need when walking the tree.
  unsigned ilevel, nremaining, k;
  double cellsize, rcelli, rcelli2, mj, rji2;
  Vector xcelli, xji, nhat;
  vector<TreeKey> newDaughters, remainingCells;
  Tree::const_iterator itr, firstLevelItr = tree.begin();

  // Walk each internal node.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {

      // State of node i.
      const double mi = mass(nodeListi, i);
      const Vector& xi = position(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& phii = potential(nodeListi, i);

      // Walk the tree.
      remainingCells = rootCell.daughters;
      while ((not remainingCells.empty())) {
        ilevel = remainingCells[0].first;
        firstLevelItr = lower_bound(tree.begin(), tree.end(), TreeKey(ilevel, CellKey(0)), TreeComparatorLessThan());
        // firstLevelItr = lower_bound(firstLevelItr, tree.end(), TreeKey(ilevel, CellKey(0)), TreeComparatorLessThan());
        nremaining = remainingCells.size();
        newDaughters = vector<TreeKey>();
        newDaughters.reserve(8*nremaining);
        cellsize = mBoxLength/(1U << ilevel);

        // Walk each of the current set of Cells.
        for (k = 0; k != nremaining; ++k) {
          itr = lower_bound(firstLevelItr, tree.end(), remainingCells[k], TreeComparatorLessThan());
          CHECK(itr != tree.end());
          const Cell& cell = itr->second;
          
          // Can we ignore this cells daughters?
          xcelli = cell.xcm - xi;
          rcelli = xcelli.magnitude();
          if (rcelli > cellsize/mOpening + cell.rcm2cc) {      // We use Barnes (1994) modified criterion.

            // Yep, treat this cells and all of its daughters as a single point.
            nhat = xcelli.unitVector();
            rcelli2 = rcelli*rcelli + soft2;
            CHECK(rcelli2 > 0.0);

            // Increment the acceleration and potential.
            DvDti += mG*cell.M/rcelli2 * nhat;
            phii -= mG*cell.M/sqrt(rcelli2);

          } else if (cell.daughters.size() == 0) {

            // This cell represents a leaf (termination of descent).  We just directly
            // add up the node properties of any nodes in the cell.
            CHECK(cell.masses.size() > 0 and
                  cell.positions.size() == cell.masses.size());
            for (unsigned k = 0; k != cell.masses.size(); ++k) {
              mj = cell.masses[k];
              xji = cell.positions[k] - xi;
              rji2 = xji.magnitude2();

              if (rji2/soft2 > 1.0e-10) {           // Screen out self-interaction.
                nhat = xji.unitVector();
                rji2 += soft2;
                CHECK(rji2 > 0.0);

                // Increment the acceleration and potential.
                DvDti += mG*mj/rji2 * nhat;
                phii -= mG*mj/sqrt(rji2);
              }
            }

          } else {

            // We need to walk further down the tree.  Add this cells daughters
            // to the next set.
            copy(cell.daughters.begin(), cell.daughters.end(), back_inserter(newDaughters));

          }
        }

        // Update the set of cells to check on the next pass.
        remainingCells = newDaughters;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Serialize a tree to a buffer of char.
//------------------------------------------------------------------------------
void
OctTreeGravity::
serialize(const OctTreeGravity::Tree& tree,
          std::vector<char>& buffer) const {
  const unsigned ncells = tree.size();
  packElement(ncells, buffer);
  for (Tree::const_iterator itr = tree.begin();
       itr != tree.end();
       ++itr) {
    packElement(itr->first.first, buffer);
    packElement(itr->first.second, buffer);
    serialize(itr->second, buffer);
  }
}

//------------------------------------------------------------------------------
// Serialize a cell to a buffer of char.
//------------------------------------------------------------------------------
void
OctTreeGravity::
serialize(const OctTreeGravity::Cell& cell,
          std::vector<char>& buffer) const {
  packElement(cell.M, buffer);
  packElement(cell.Mglobal, buffer);
  packElement(cell.xcm, buffer);
  packElement(cell.rcm2cc, buffer);
  packElement(cell.daughters, buffer);
  packElement(cell.masses, buffer);
  packElement(cell.positions, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a tree from a buffer of char.
//------------------------------------------------------------------------------
void
OctTreeGravity::
deserialize(OctTreeGravity::Tree& tree,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unsigned ncells;
  TreeKey key;
  Cell cell;
  unpackElement(ncells, bufItr, endItr);
  for (unsigned i = 0; i != ncells; ++i) {
    cell.daughters = vector<TreeKey>();
    cell.masses = vector<double>();
    cell.positions = vector<Vector>();
    unpackElement(key.first, bufItr, endItr);
    unpackElement(key.second, bufItr, endItr);
    deserialize(cell, bufItr, endItr);
    tree.push_back(make_pair(key, cell));
  }
}

//------------------------------------------------------------------------------
// Deserialize a cell from a buffer of char.
//------------------------------------------------------------------------------
void
OctTreeGravity::
deserialize(OctTreeGravity::Cell& cell,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unpackElement(cell.M, bufItr, endItr);
  unpackElement(cell.Mglobal, bufItr, endItr);
  unpackElement(cell.xcm, bufItr, endItr);
  unpackElement(cell.rcm2cc, bufItr, endItr);
  unpackElement(cell.daughters, bufItr, endItr);
  unpackElement(cell.masses, bufItr, endItr);
  unpackElement(cell.positions, bufItr, endItr);
}

//------------------------------------------------------------------------------
// Define our static members.
//------------------------------------------------------------------------------
unsigned OctTreeGravity::num1dbits = 21U;
uint64_t OctTreeGravity::max1dKey = 1U << OctTreeGravity::num1dbits;
uint64_t OctTreeGravity::xkeymask = (1U << OctTreeGravity::num1dbits) - 1U;
uint64_t OctTreeGravity::ykeymask = OctTreeGravity::xkeymask << OctTreeGravity::num1dbits;
uint64_t OctTreeGravity::zkeymask = OctTreeGravity::ykeymask << OctTreeGravity::num1dbits;

}
}

#include "OctTreeGravityInline.hh"
