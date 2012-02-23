//---------------------------------Spheral++----------------------------------//
// NBodyGravity -- A simple implementation of the N-body gravity algorithm.
//
//! \author $Author: mikeowen $
//! \version $Revision: 3174 $
//! \date $Date: 2009-07-22 14:23:24 -0700 (Wed, 22 Jul 2009) $
//
//----------------------------------------------------------------------------//
#ifndef NBODYGRAVITY_HH
#define NBODYGRAVITY_HH

#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace GravitySpace {

template <typename Dimension>
class NBodyGravity: public PhysicsSpace::GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;

  //! Constructor.
  //! \param plummerSofteningLength -- The Plummer Softening Length for the model.
  //! \param maxDeltaVelocity -- Maximum factor by which the velocity can be changed by an acceleration per timestep.
  NBodyGravity(const double plummerSofteningLength,
               const double maxDeltaVelocity,
               const double G);

  //! Destructor.
  virtual ~NBodyGravity();

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

  //! Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The current softening length.
  double softeningLength() const;
  void softeningLength(const double x);

private:
  
  //! The gravitational potential of the particles.
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;

  //! The total potential energy of the particles.  Also mutable.
  mutable Scalar mExtraEnergy;

  //! The maximum acceleration magnitude during the last time step.
  mutable Scalar mOldMaxAcceleration;
  
  //! The maximum velocity magnitude during the last time step.
  mutable Scalar mOldMaxVelocity;

  //! The maximum allowed change in velocity, by factors of velocity.
  Scalar mMaxDeltaVelocityFactor;

  //! The Plummer softening length.
  Scalar mSofteningLength;
  
  //! The gravitational constant.
  Scalar mG;
  
  // Default constructor -- disabled.
  NBodyGravity();

  // Copy constructor -- disabled.
  NBodyGravity(const NBodyGravity&);

  // Assignment operator -- disabled.
  NBodyGravity& operator=(const NBodyGravity&);

}; // end class NBodyGravity

} // end namespace GravitySpace

} // end namespace Spheral

#endif
