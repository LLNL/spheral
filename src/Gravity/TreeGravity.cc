//---------------------------------Spheral++----------------------------------//
// TreeGravity implementation.
//
//! \author $Author: mikeowen $
//! \version $Revision: 4012 $
//! \date $Date: 2011-02-20 10:58:32 -0800 (Sun, 20 Feb 2011) $
//----------------------------------------------------------------------------//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "TreeGravity.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "DBC.hh"
#include "Material/PhysicalConstants.hh"

#include "libs.hh"
#include "classes.hh"
#include "headers.hh"

namespace Spheral {
namespace GravitySpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
TreeGravity::
TreeGravity(const double G,
            const double openeing):
  mG(G),
  mOpening(opening),
  mPotential(FieldList<Dim<3>, Scalar>::Copy),
  mExtraEnergy(0.0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
TreeGravity::
~TreeGravity() {
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
void 
TreeGravity::
evaluateDerivatives(const Dim<3>::Scalar time,
                    const Dim<3>::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dim<3> >& state,
                    StateDerivatives<Dim<3> >& derivs) const {
  using namespace NodeSpace;

  // Access to pertinent fields in the database.
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

  // We scale the positions to be in the unit cube, and the masses such that G
  // is one inside Tree.
  const double boxlength = mXmax.x() - mXmin.x();
  VERIFY(fuzzyEqual(mXmax.y() - mXmin.y(), boxlength, 1.0e-10));
  VERIFY(fuzzyEqual(mXmax.z() - mXmin.z(), boxlength, 1.0e-10));
  VERIFY(boxlength > 0.0);
  const double lscale = 1.0/boxlength;
  const double mscale = lscale*lscale*lscale/mG;

  // Create the Tree memory struct, and fill in some of it's parameters.
  // For now we will scale the input to Tree for unit length and unit total mass.
  TreeSpace::Tree_Memory* pmemory = generateTreeMemory(numNodes,
                                                                mPeriodic,
                                                                mNgrid,
                                                                mMinHighParticles,
                                                                mPadding,
                                                                1.0);

  // Create the Tree class.
  TreeSpace::Tree* pfrac = new TreeSpace::Tree(*pmemory);
  pmemory->p_fractal = pfrac;
  pfrac->particle_list.resize(pmemory->number_particles);

  // Create the memory for Tree's particles.
  vector<TreeSpace::Particle> particles(numNodes);
  for (unsigned i = 0; i != numNodes; ++i) {
    particles[i].space_resize(3);    // 3 => pos(x, y, z)
    particles[i].field_resize(4);    // 4 => fields(pot, fx, fy, fz)
  }

  // Copy Spheral's particle information to Tree's particles.
  unsigned j = 0;
  double mtot = 0.0;
  vector<double> ppos(3);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {
      CHECK(j < numNodes);
      pfrac->particle_list[j] = &particles[j];
      const Vector xi = (position(nodeListi, i) - mXmin)/lscale;
      ppos[0] = max(0.0, min(1.0, xi.x()));
      ppos[1] = max(0.0, min(1.0, xi.y()));
      ppos[2] = max(0.0, min(1.0, xi.z()));
      particles[j].set_mass(mass(nodeListi, i)/mscale);
      particles[j].set_pos(ppos);
      mtot += particles[j].get_mass();
      ++j;
    }
  }
  pmemory->total_mass = mtot;

  // Invoke the gravity solver.
  pfrac->timing(-2,0);
  pfrac->timing(-1,29);
  TreeSpace::fractal_gravity(*pfrac, *pmemory);
  pfrac->timing(1,29);
  pfrac->timing(0,0);

  // Read the result back to Spheral's data structures.
  vector<double> f(4);
  j = 0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {
      CHECK(j < numNodes);
      DxDt(nodeListi, i) += velocity(nodeListi, i);

      // Extract field (pot, fx, fy, fz) from the Tree particle.
      particles[j].get_field_pf(f);

      // Update the potential energy.
      mPotential(nodeListi, i) = f[0] * mscale*lscale*lscale;
      mExtraEnergy += mass(nodeListi, i) * mPotential(nodeListi, i);

      // Update the acceleration.
      DvDt(nodeListi, i) += Vector(f[1] * lscale,
                                   f[2] * lscale,
                                   f[3] * lscale);

      // Capture the maximum acceleration and velocity magnitudes.
      const Scalar accelMagnitude = DvDt(nodeListi, i).magnitude();
      if (mOldMaxAcceleration < accelMagnitude) {
        mOldMaxAcceleration = accelMagnitude + 1e-10;
        CHECK(mOldMaxAcceleration != 0.0);
      }

      const Scalar velocityMagnitude = velocity(nodeListi, i).magnitude();
      if (mOldMaxVelocity < velocityMagnitude) {
        mOldMaxVelocity = velocityMagnitude;
      }
    }
  }

  // Clean up.
  delete pfrac;
  delete pmemory;
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
void 
TreeGravity::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}

//------------------------------------------------------------------------------
// Vote on a time step.  We should fill in a sqrt(G/rho) type thing here!
//------------------------------------------------------------------------------
TreeGravity::TimeStepType
TreeGravity::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

  // The maximum change in our velocity for the next time cycle is given 
  // by the mMaxDeltaVelocityFactor plus the max velocity and the max 
  // accelation from the last cycle.
  const double deltat = (mOldMaxVelocity*mOldMaxAcceleration/(mOldMaxAcceleration*mOldMaxAcceleration + 1.0e-10)) * mMaxDeltaVelocityFactor;

  stringstream reasonStream;
  reasonStream << "velocity: " << mOldMaxVelocity
               << ", acceleration: " << mOldMaxAcceleration
               << "dt = f*v/a: " << deltat
               << ends;
  return TimeStepType(deltat, reasonStream.str());
}

//------------------------------------------------------------------------------
// extraEnergy
//------------------------------------------------------------------------------
TreeGravity::Scalar 
TreeGravity::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
const FieldList<Dim<3>, TreeGravity::Scalar>&
TreeGravity::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
double
TreeGravity::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// opening
//------------------------------------------------------------------------------
double
TreeGravity::
opening() const {
  return mOpening;
}

}
}
