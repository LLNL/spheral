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
    generateFractalMemory(const unsigned numNodes,
			  const bool periodic,
			  const unsigned ngrid,
			  const unsigned minHighParticles,
			  const unsigned padding,
			  const double boxlength) {
      using namespace FractalSpace;
      
      Fractal_Memory* pmemory = new FractalSpace::Fractal_Memory;
      pmemory->periodic = periodic;
      pmemory->grid_length = ngrid;
      pmemory->minimum_number = minHighParticles;
      pmemory->padding = padding;
      pmemory->box_length =boxlength;
      pmemory->number_particles = numNodes;
      
      pmemory->debug=true;
      pmemory->new_points_gen=9;
      pmemory->remember_points=false;
      pmemory->number_steps_total=903;
      pmemory->redshift_start=99.0;
      pmemory->max_particles=300000;
      pmemory->omega_0=0.3;
      pmemory->omega_lambda=0.7;
      pmemory->h=0.7;
      pmemory->steps=-1;
      pmemory->random_gen=54321;
      // sets your values for parameters
      pmemory->amnesia=true; // (true) forget everything after you are done. (false) remember everything.
      pmemory->mind_wipe=false; // (true) delete everything and then come back without calculating anything.
      pmemory->fixed_potential=false; // (true) use the fixed potential.
      pmemory->calc_shear=false;// (true) if we calculate shear of force field
      pmemory->calc_density_particle=false;
      pmemory->do_vel=false;
      pmemory->start_up=false;
      pmemory->halo=false;
      pmemory->halo_fixed=pmemory->halo_fixed && !pmemory->periodic;
      pmemory->length_ratio=1;
      pmemory->omega_start=Omega(pmemory->omega_0,pmemory->omega_lambda,pmemory->redshift_start);
      pmemory->lambda_start=Lambda(pmemory->omega_0,pmemory->omega_lambda,pmemory->redshift_start);
      pmemory->sigma_initial=pmemory->sigma_0*Growth(pmemory->omega_0,pmemory->omega_lambda,pmemory->redshift_start);
      pmemory->time=Age_of_the_universe(pmemory->omega_start,pmemory->lambda_start,0.0);
      pmemory->total_mass=1.0;
      //
      pmemory->crash_levels=5;
      pmemory->crash_pow=2.0;
      pmemory->density_crash=5.5;
      pmemory->splits=0;
      //
      pmemory->masks=0;
      //
      pmemory->masks_init=0;

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
      const double lunscale = 1.0/lscale;
      const double munscale = 1.0/mscale;

      // Create the Fractal memory struct, and fill in some of it's parameters.
      // For now we will scale the input to Fractal for unit length and unit total mass.
      FractalSpace::Fractal_Memory* pmemory = generateFractalMemory(numNodes,
								    mPeriodic,
								    mNgrid,
								    mMinHighParticles,
								    mPadding,
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
      delete pmemory;
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
