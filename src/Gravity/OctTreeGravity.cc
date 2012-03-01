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
#include "Utilities/boundingBox.hh"
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
               const double opening,
               const double softeningLength,
               const double ftimestep):
  mG(G),
  mOpening(opening),
  mSofteningLength(softeningLength),
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
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
void 
OctTreeGravity::
evaluateDerivatives(const Dim<3>::Scalar time,
                    const Dim<3>::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dim<3> >& state,
                    StateDerivatives<Dim<3> >& derivs) const {

  const double soft2 = mSofteningLength*mSofteningLength;

  // Access the pertinent fields in the database.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  const unsigned numNodes = mass.numInternalNodes();

  // Get the acceleration and position change vectors we'll be modifying.
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;

  // We'll always be starting with the daughters of the root level.
  const unsigned numLevels = mTree.size();
  CHECK(numLevels >= 1);
  const Cell& rootCell = mTree[0].begin()->second;

  // Walk each internal node.
  TreeLevel::const_iterator cellItr;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {

      // State of node i.
      const double mi = mass(nodeListi, i);
      const Vector& xi = position(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& phii = mPotential(nodeListi, i);

      // Set the position update to Lagrangian.
      DxDt(nodeListi, i) += velocity(nodeListi, i);

      // Walk the tree.
      unsigned ilevel = 0;
      vector<CellKey> remainingCells = rootCell.daughters;
      while ((not remainingCells.empty()) and ++ilevel < numLevels) {
        vector<CellKey> newDaughters;
        const double cellsize = mBoxLength/(1U << ilevel);

        // Walk each of the current set of Cells.
        for (vector<CellKey>::const_iterator keyItr = remainingCells.begin();
             keyItr != remainingCells.end();
             ++keyItr) {
          cellItr = mTree[ilevel].find(*keyItr);
          CHECK(cellItr != mTree[ilevel].end());
          const Cell& cell = cellItr->second;
          
          // Can we ignore this cells daughters?
          const Vector xcelli = cell.xcm - xi;
          const double rcelli = xcelli.magnitude();
          if (rcelli > cellsize/mOpening + cell.rcm2cc) {      // We use Barnes (1994) modified criterion.

            // Yep, treat this cells and all of its daughters as a single point.
            const Vector nhat = xcelli.unitVector();
            const double rcelli2 = rcelli*rcelli + soft2;
            CHECK(rcelli2 > 0.0);

            // Increment the acceleration and potential.
            DvDti += mG*cell.M/rcelli2 * nhat;
            phii -= mG*cell.M/sqrt(rcelli2);

//             cerr << "node "<< i << " interacting with cell (" << ilevel << " " << *keyItr << ") : " << cell.M << " " << rcelli2 << endl;

          } else if (cell.daughters.size() == 0) {

            // This cell represents a leaf (termination of descent.  We just directly
            // add up the node properties of any nodes in the cell.
            CHECK(cell.members.size() > 0);
            for (unsigned k = 0; k != cell.members.size(); ++k) {
              const unsigned nodeListj = cell.members[k].first;
              const unsigned j = cell.members[k].second;

              if (nodeListj != nodeListi or j != i) {           // Screen out self-interaction.
                const Vector xji = position(nodeListj, j) - xi;
                const Vector nhat = xji.unitVector();
                const double rji2 = xji.magnitude2() + soft2;
                CHECK(rji2 > 0.0);

                // Increment the acceleration and potential.
                const double mj = mass(nodeListj, j);
                DvDti += mG*mj/rji2 * nhat;
                phii -= mG*mj/sqrt(rji2);

//                 cerr << "node " << i << " interacting with point (" << ilevel << " " << *keyItr << ") : (" << nodeListj << " " << j << ") : " << mj << " " << rji2 << endl;
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
  boundingBox(position, mXmin, mXmax, false, false);
  CHECK(mXmin.x() <= mXmax.x());
  CHECK(mXmin.y() <= mXmax.y());
  CHECK(mXmin.z() <= mXmax.z());
  mBoxLength = (mXmax - mXmin).maxAbsElement();
  CHECK(mBoxLength > 0.0);

  // Walk all the internal nodes and add them to the tree.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const size_t n = mass[nodeListi]->numInternalElements();
    for (size_t i = 0; i != n; ++i) {
      this->addNodeToTree(nodeListi, i, mass(nodeListi, i), position(nodeListi, i));
    }
  }

  // Make a final pass over the cells and fill in the distance between the 
  // center of mass and the geometric center.
  // We also squirrel away the maximum effective cell density for our timestep
  // determination.
  CellKey ckey, ix, iy, iz;
  double cellsize, cellvol;
  mMaxCellDensity = 0.0;
  for (unsigned ilevel = 0; ilevel != mTree.size(); ++ilevel) {
    cellsize = mBoxLength/(1U << ilevel);
    cellvol = cellsize*cellsize*cellsize;
    for (TreeLevel::iterator itr = mTree[ilevel].begin();
         itr != mTree[ilevel].end();
         ++itr) {
      ckey = itr->first;
      Cell& cell = itr->second;
      extractCellIndices(ckey, ix, iy, iz);

      // Update the distance between the cell's center of mass and geometric center.
      cell.rcm2cc = (cell.xcm - (mXmin + Vector((ix + 0.5)*cellsize,
                                                (iy + 0.5)*cellsize,
                                                (iz + 0.5)*cellsize))).magnitude();
      CHECK(cell.rcm2cc < 1.74*cellsize);

      // Update the maximum effective cell density.
      mMaxCellDensity = max(mMaxCellDensity, cell.M/cellvol);
    }
  }
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
  const double dt = mftimestep * sqrt(mG/mMaxCellDensity);

  stringstream reasonStream;
  reasonStream << "OctTreeGravity: sqrt(/(G rho)) = sqrt(1/("
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
  ss << "Tree : nlevels = " << mTree.size() << "\n";
  for (unsigned ilevel = 0; ilevel != mTree.size(); ++ilevel) {
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << mTree[ilevel].size() << "\n";
    for (TreeLevel::const_iterator itr = mTree[ilevel].begin();
         itr != mTree[ilevel].end();
         ++itr) {
      key = itr->first;
      const Cell& cell = itr->second;
      extractCellIndices(key, ix, iy, iz);
      ss << "    Cell key=" << key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << "\n"
         << "         xcm=" << cell.xcm << " rcm2cc=" << cell.rcm2cc << " M=" << cell.M << "\n"
         << "         daughters = ( ";
      for (vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) ss << *ditr << " ";
      ss << ")\n"
         << "         nodes = [";
      for (vector<NodeID>::const_iterator nitr = cell.members.begin();
           nitr != cell.members.end();
           ++nitr) ss << " (" << nitr->first << " " << nitr->second <<")";
      ss <<" ]\n";
    }
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
