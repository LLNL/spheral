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
#include "classes.hh"

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
  //! \param ngrid -- length of the fundamental grid.
  //! \param maxDeltaVelocity -- Maximum factor by which the velocity can be changed by an 
  //!                            acceleration per timestep.
  FractalGravity(const double G,
                 const Vector& xmin,
                 const Vector& xmax,
                 const bool periodic,
                 const unsigned gridLength,
                 const double maxDeltaVelocity);

  //! Destructor.
  virtual ~FractalGravity();

  //! This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  //! Initializations on problem start up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& db) override;

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const override;

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
  unsigned gridLength() const;

private:
  //! Gravitational constant
  double mG;

  //! Box boundaries
  Vector mXmin, mXmax;

  //! Choose periodic or isolated boundaries
  bool mPeriodic;

  //! Length of the grid for the FFT box
  unsigned mGridLength;
  
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

  //! The FractalGravity memory blob
  FractalSpace::Fractal_Memory* mFractalMemoryPtr;

  //! Special communicator for Fractal.
  MPI_Comm mFractalComm;

  // Disabled methods.
  FractalGravity();
  FractalGravity(const FractalGravity&);
  FractalGravity& operator=(const FractalGravity&);
};

}
}

#endif
