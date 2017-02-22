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
#include "Utilities/DBC.hh"
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
// This is an internal convenience method to create Jens' Fractal_Memory struct.
// I have no idea what most of the parameters in this thing should be, so I'm 
// copying one of Jens' examples (fractal_memory_parameters_gal).  
// This should be reviewed by Jens!
//------------------------------------------------------------------------------
FractalSpace::Fractal_Memory* 
generateFractalMemory(const int numNodes,
		      const bool periodic,
		      const unsigned ngrid,
		      const unsigned nlevelmax,
		      const unsigned minHighParticles,
		      const unsigned padding,
		      const double tolHypre,
		      const unsigned maxitsHypre,
		      const bool fractalDebug,
		      const unsigned FractalNodes0,
		      const unsigned FractalNodes1,
		      const unsigned FractalNodes2,
		      const string BaseDirectory,
		      const string RunIdentifier){
  using namespace FractalSpace;

  Fractal_Memory* pmemory = new FractalSpace::Fractal_Memory;
  pmemory->number_particles = numNodes;
  pmemory->periodic = periodic;
  pmemory->grid_length = ngrid;
  pmemory->level_max=nlevelmax;
  pmemory->minimum_number = minHighParticles;
  pmemory->padding = padding;
  pmemory->epsilon_sor=tolHypre;
  pmemory->maxits=maxitsHypre;
  pmemory->debug=fractaldebug;
  pmemory->FractalNodes0=FractalNodes0;
  pmemory->FractalNodes1=FractalNodes1;
  pmemory->FractalNodes1=FractalNodes2;
  pmemory->FractalNodes=FractalNodes0*FractalNodes1*FractalNodes2;
  pmemory->MPIrun=pmemory.MPIrun || pmemory->FractalNodes > 1;
  pmemory->BaseDirectory=BaseDirectory;
  pmemory->RUN=RunIdentifier;

  // That's it.
  return pmemory;
}

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
	       const double epsilon_poisson,
	       const unsigned maxits_poisson,
	       const bool fractal_debug,
	       const unsigned FractalNodes0,
	       const unsigned FractalNodes1,
	       const unsigned FractalNodes2,
	       const double maxDeltaVelocity);
  mG(G),
  mXmin(xmin),
  mXmax(xmax),
  mPeriodic(periodic),
  mNgrid(ngrid),
  mNlevelmax(nlevelmax),
  mMinHighParticles(minHighParticles),
  mPadding(padding),
  mEpsilonPoisson(epsilon_poisson),
  mMaxitsPoisson(maxits_poisson),
  mFractalDebug(fractal_debug),
  mFractalNodes0(FractalNodes0),
  mFractalNodes1(FractalNodes1),
  mFractalNodes2(FractalNodes2),
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
  const double lunscale = 1.0/lscale;
  const double munscale = 1.0/mscale;
  
  // Create the Fractal memory struct, and fill in some of it's parameters.
  // For now we will scale the input to Fractal for unit length and unit total mass.
  FractalSpace::Fractal_Memory* pmemory = generateFractalMemory(
								numNodes,
								mPeriodic,
								Mngrid,
								nlevelmax,
								minHighParticles,
								padding,
								tolHypre,
								maxitsHypre,
								fractalDebug,
								FractalNodes0,
								FractalNodes1,
								FractalNodes2,
								BaseDirectory,
								RunIdentifier)
    bool MR=pmemory->MPIrun;
  int GR=pmemory->grid_length;
  bool PR=pmemory->periodic;
  int NP=pmemory->number_particles;
  Mess* p_mess=new Mess(MR,GR,PR,NP);
  pmemory->p_mess=p_mess;
  string BD=pmemory->BaseDirectory;
  int FR=p_mess->FractalRank;
  string RUN=pmemory->RUN;
  File* p_file=new File(BD,FR,RUN);
  pmemory->p_file=p_file;
  fractal_force_init(pmemory);
                                                                1.0);

  // Create the Fractal class.
  FractalSpace::Fractal* pfrac = new FractalSpace::Fractal(*pmemory);
  pfrac->set_number_particles(numNodes);
  pmemory->p_fractal = pfrac;
  pfrac->particle_list.resize(pmemory->number_particles);

  // Create the memory for Fractal's particles.
  vector<FractalSpace::Particle> particles(numNodes);
  for (unsigned i = 0; i != numNodes; ++i) {
    particles[i].space_resize(3);    // 3 => pos(x, y, z)
    particles[i].field_resize(4);    // 4 => fields(pot, fx, fy, fz)
  }

  // Copy Spheral's particle information to Fractal's particles.
  unsigned j = 0;
  double mtot = 0.0;
  vector<double> ppos(3);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {
      CHECK(j < numNodes);
      pfrac->particle_list[j] = &particles[j];
      const Vector xi = (position(nodeListi, i) - mXmin)*lscale;
      ppos[0] = max(0.0, min(1.0, xi.x()));
      ppos[1] = max(0.0, min(1.0, xi.y()));
      ppos[2] = max(0.0, min(1.0, xi.z()));
      particles[j].set_mass(mass(nodeListi, i)*mscale);
      particles[j].set_pos(ppos);
      mtot += particles[j].get_mass();
      ++j;
    }
  }
  pmemory->total_mass = mtot;

  // Invoke the gravity solver.
  pfrac->timing(-2,0);
  pfrac->timing(-1,29);
  FractalSpace::fractal_gravity(*pfrac, *pmemory);
  pfrac->timing(1,29);
  pfrac->timing(0,0);

  // Read the result back to Spheral's data structures.
  vector<double> f(4);
  j = 0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {
      CHECK(j < numNodes);
      DxDt(nodeListi, i) += velocity(nodeListi, i);

      // Extract field (pot, fx, fy, fz) from the Fractal particle.
      particles[j].get_field_pf(f);

      // Update the potential energy.
      mPotential(nodeListi, i) = f[0] * munscale*lunscale*lunscale;
      mExtraEnergy += mass(nodeListi, i) * mPotential(nodeListi, i);

      // Update the acceleration.
      DvDt(nodeListi, i) += Vector(f[1] * lunscale,
                                   f[2] * lunscale,
                                   f[3] * lunscale);

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
  //  delete pmemory;
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
void 
FractalGravity::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}

//------------------------------------------------------------------------------
// Vote on a time step.  We should fill in a sqrt(G/rho) type thing here!
//------------------------------------------------------------------------------
FractalGravity::TimeStepType
FractalGravity::
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
FractalGravity::Scalar 
FractalGravity::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
const FieldList<Dim<3>, FractalGravity::Scalar>&
FractalGravity::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
double
FractalGravity::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
FractalGravity::Vector
FractalGravity::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
FractalGravity::Vector
FractalGravity::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// periodic
//------------------------------------------------------------------------------
bool
FractalGravity::
periodic() const {
  return mPeriodic;
}

//------------------------------------------------------------------------------
// ngrid
//------------------------------------------------------------------------------
unsigned
FractalGravity::
ngrid() const {
  return mNgrid;
}

//------------------------------------------------------------------------------
// nlevelmax
//------------------------------------------------------------------------------
unsigned
FractalGravity::
nlevelmax() const {
  return mNlevelmax;
}

//------------------------------------------------------------------------------
// minHighParticles
//------------------------------------------------------------------------------
unsigned
FractalGravity::
minHighParticles() const {
  return mMinHighParticles;
}

//------------------------------------------------------------------------------
// padding
//------------------------------------------------------------------------------
unsigned
FractalGravity::
padding() const {
  return mPadding;
}

}
}
