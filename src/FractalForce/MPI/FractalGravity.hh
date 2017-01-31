//---------------------------------Spheral++----------------------------------//
// FractalGravity.
//
// A Spheral++ interface to Jens Villumsen's fractal gravity HPM Poisson solver.
//
//! \author $Author: mikeowen $
//! \version $Revision: 3174 $
//! \date $Date: 2009-07-22 14:23:24 -0700 (Wed, 22 Jul 2009) $
//
//----------------------------------------------------------------------------//
#ifndef __Spheral_FractalGravity_interface__
#define __Spheral_FractalGravity_interface__

#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace GravitySpace {

class FractalGravity: public PhysicsSpace::GenericBodyForce<Dim<3> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;

  //! Constructor.
  //! \param G -- the gravitational constant.
  //! \param xmin -- lower left corner of computational cube.
  //! \param xmax -- upper right corner of computational cube.
  //! \param periodic -- should the problem be solved in periodic (true) 
  //!                    or isolated (false) geometry.
  //! \param ngrid -- length of the fundamental grid. Must be even.
  //! \param nlevelmax -- maximum number of levels for refinement.
  //! \param minHighParticles -- minimum number of particles belonging to a high point.
  //! \param padding --(-1) => resolution jumps by factors of 2 (cheap)
  //!                   (1) => each high point fully padded (expensive)
  //!                   (0) => no padding (not recommended)
  //! \param tolHypre -- accuracy of hypre solver (1.0e-7 recommended)
  //! \param maxitsHypre -- max number of iterations in Hypre (20 recommended)
  //! \param fractalDebug -- debug run (true) extra IO
  //! \param FractalNodes0 -- Number of processors in x-direction
  //! \param FractalNodes1 -- Number of processors in y-direction
  //! \param FractalNodes2 -- Number of processors in z-direction
  //! \param BaseDirectory -- (string) Base directory for Fractal output
  //! \param RunIdentifier -- (string) Identifier
  //! \param maxDeltaVelocity -- Maximum factor by which the velocity can be changed by an 
  //!                            acceleration per timestep.
  FractalGravity(const double G,
                 const Vector& xmin,
                 const Vector& xmax,
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
		 const string RunIdentifier,
                 const double maxDeltaVelocity);

  //! Destructor.
  virtual ~FractalGravity();

  //! This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Initializations on problem start up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& db);

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  //! Return the gravitational potential created by the particle distribution.
  const FieldSpace::FieldList<Dimension, Scalar>& potential() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The lower left corner of the computational cube.
  Vector xmin() const;

  //! The upper right corner of the computational cube.
  Vector xmax() const;

  //! Should the box be treated as periodic or isolated?
  bool periodic() const;

  //! Number of grid cells per side on the fundamental level.
  unsigned ngrid() const;

  //! Maximum number of levels of refinement.
  unsigned nlevelmax() const;

  //! Minimum number of particles to form a high point.
  unsigned minHighParticles() const;

  //! padding flag:
  //!    -1 => resolution jumps by factors of 2 (cheap)
  //!     0 => no padding (not recommended)
  //!     1 => each high point fully padded (expensive)
  unsigned padding() const;

  //! accuracy in Hypre Poisson Solver (1.0e-7 standard)
  double tolHypre() const;

  //! iterations in Hypre Poisson Solver (20 standard)
  unsigned maxitsHypre() const;

  //! debug run true/false
  bool fractalDebug() const;

  //! Number of nodes in x-direction
  unsigned FractalNodes0() const;

  //! Number of nodes in y-direction
  unsigned FractalNodes1() const;

  //! Number of nodes in z-direction
  unsigned FractalNodes2() const;

  //! Basedirectory for Fractal output, I use "/p/lscratchc/jensv"
  string BaseDirectory() const;

  //! Identifier for run
  string RunIdentifier();
private:
  double mG;
  Vector mXmin, mXmax;
  bool mPeriodic;
  unsigned mNgrid, mNlevelmax, mMinHighParticles, mPadding;
  
  //! The maximum allowed change in velocity, by factors of velocity.
  Scalar mMaxDeltaVelocityFactor;

  //! The gravitational potential of the particles.
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;

  //! The total potential energy of the particles.  Also mutable.
  mutable Scalar mExtraEnergy;

  //! The maximum acceleration magnitude during the last time step.
  mutable Scalar mOldMaxAcceleration;
  
  //! The maximum velocity magnitude during the last time step.
  mutable Scalar mOldMaxVelocity;

  // Default constructor -- disabled.
  FractalGravity();

  // Copy constructor -- disabled.
  FractalGravity(const FractalGravity&);

  // Assignment operator -- disabled.
  FractalGravity& operator=(const FractalGravity&);

};

}
}

#endif
