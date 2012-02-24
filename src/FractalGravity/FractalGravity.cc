//---------------------------------Spheral++----------------------------------//
// FractalGravity implementation.
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

#include "FractalGravity.hh"
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

namespace Spheral {
namespace GravitySpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
FractalGravity::
FractalGravity(const double G,
               const Vector& xmin,
               const Vector& xmax,
               const bool periodic,
               const unsigned ngrid,
               const unsigned nlevelmax,
               const unsigned minHighParticles,
               const unsigned padding,
               const double maxDeltaVelocity):
  mG(G),
  mXmin(xmin),
  mXmax(xmax),
  mPeriodic(periodic),
  mNgrid(ngrid),
  mNlevelmax(nlevelmax),
  mMinHighParticles(minHighParticles),
  mPadding(padding),
  mMaxDeltaVelocityFactor(maxDeltaVelocity),
  mPotential(FieldList<Dim<3>, Scalar>::Copy),
  mExtraEnergy(0.0),
  mOldMaxAcceleration(0.0),
  mOldMaxVelocity(0.0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
FractalGravity::
~FractalGravity() {
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
void 
FractalGravity::
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
  // is one inside Fractal.
  const double boxlength = mXmax.x() - mXmin.x();
  VERIFY(fuzzyEqual(mXmax.y() - mXmin.y(), boxlength, 1.0e-10));
  VERIFY(fuzzyEqual(mXmax.z() - mXmin.z(), boxlength, 1.0e-10));
  VERIFY(boxlength > 0.0);
  const double lscale = 1.0/boxlength;
  const double mscale = lscale*lscale*lscale/mG;

  // Create the Fractal memory.
  FractalSpace::Fractal_memory* pmemory = new Fractal_Memory;

  // Fill in the Fractal_memory struct.
  // For now we will scale the input to Fractal for unit length and unit total mass.
  pmemory->periodic = mPeriodic;
  pmemory->grid_length = mNgrid;
  pmemory->minimum_number = mMinHighParticles;
  pmemory->padding = mPadding;
  pmemory->box_length = 1.0;
  pmemory->number_particles = numNodes;

  // Create the Fractal class.
  FractalSpace::Fractal* pfrac = new Fractal(pmemory);
  pmemory->p_fractal = pfrac;
  CHECK(pfrac->number_particles == numNodes);

  // Walk the particles and fill in Fractal's particle structures.
  vector<FractalSpace::Particle> particles;
  particles.resize(numNodes);
  unsigned j = 0;
  double mtot = 0.0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi].numInternalElements(); ++i) {
      CHECK(j < numNodes);
      pfrac->particle_list[j] = &particles[j];
      const Vector xi = (position(nodeListi, i) - mXmin)/lscale;
      particles[j].mass = mass(nodeListi, i)/mscale;
      particles[j].phase_space[0] = max(0.0, min(1.0, xi.x()));
      particles[j].phase_space[1] = max(0.0, min(1.0, xi.y()));
      particles[j].phase_space[2] = max(0.0, min(1.0, xi.z()));
      mtot += particles[j].mass;
      ++j;
    }
  }
  pmemory->total_mass = mtot;

  // Invoke the gravity solver.
  fractal_gravity(pfrac, pmemory);

  // Read the result back to Spheral's data structures.
  unsigned ifield = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator iitr = dataBase.fluidNodeListBegin();
       iitr != dataBase.fluidNodeListEnd();
       ++iitr, ++ifield) {
    const NodeList<Dimension>& nodeListi = **iitr;
    const unsigned firstGhostNodei = nodeListi.firstGhostNode();

    for (unsigned i = 0; i != firstGhostNodei; ++i) {
      unsigned jfield = 0;

      // Set the position derivative.
      DxDt(ifield, i) += velocity(ifield, i);

      // Get a reference to the acceleration vector.
      Vector& acceleration = DvDt(ifield, i);
    
      // Zero out the potential of this particle.
      mPotential(ifield, i) = 0.0;

      for (typename DataBase<Dimension>::ConstFluidNodeListIterator jitr = dataBase.fluidNodeListBegin();
           jitr != dataBase.fluidNodeListEnd();
           ++jitr, ++jfield) {
        const NodeList<Dimension>& nodeListj = **jitr;
        const unsigned nj = nodeListj.numNodes();
        for (unsigned j = 0; j != nj; ++j) {

          // Particles can't self-interact, silly!
          if (ifield != jfield or i != j) {

            // Contribute to the acceleration of this particle.
            const Vector r = position(ifield, i) - position(jfield, j);

            CHECK(r.magnitude2() != 0.0);
            Vector rHat = r.unitVector();
            Scalar distance2 = r.magnitude2() + softeningLength2;
            CHECK(distance2 != 0.0);
            acceleration -= mG * mass(jfield, j) * rHat / distance2;

            // Also sum up contributions to the potential and 
            // total potential energy.
            mPotential(ifield, i) += mG * mass(jfield, j) / std::sqrt(distance2);
          }
        }
      }

      mExtraEnergy += mG * mass(ifield, i) * mPotential(ifield, i);

      // Capture the maximum acceleration and velocity magnitudes.
      const Scalar accelMagnitude = acceleration.magnitude();
      if (mOldMaxAcceleration < accelMagnitude) {
        mOldMaxAcceleration = accelMagnitude + 1e-10;
        CHECK(mOldMaxAcceleration != 0.0);
      } // end if

      const Scalar velocityMagnitude = velocity(ifield, i).magnitude();
      if (mOldMaxVelocity < velocityMagnitude) {
        mOldMaxVelocity = velocityMagnitude;
      }
    }
  }
}

//------------------------------------------------------------------------------
// initialize()
//------------------------------------------------------------------------------
template <typename Dimension>
void 
FractalGravity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename FractalGravity<Dimension>::TimeStepType
FractalGravity<Dimension>::
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

//------------------------------------------------------------------------------
template <typename Dimension>
typename FractalGravity<Dimension>::Scalar 
FractalGravity<Dimension>::
extraEnergy() const {
  return mExtraEnergy;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
const FieldList<Dimension, typename FractalGravity<Dimension>::Scalar>&
FractalGravity<Dimension>::
potential() const {
  return mPotential;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
bool 
FractalGravity<Dimension>::
valid() const {
  // We're always valid, sir!  (This is crap, but we can make no other 
  // assumptions right now.)
  return true;
}

//------------------------------------------------------------------------------
template <typename Dimension>
double
FractalGravity<Dimension>::
G() const {
  return mG;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
double
FractalGravity<Dimension>::
softeningLength() const {
  return mSofteningLength;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
void
FractalGravity<Dimension>::
softeningLength(const double x) {
  VERIFY(x >= 0.0);
  mSofteningLength = x;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
}
}
