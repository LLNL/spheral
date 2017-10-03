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

#include <vector>

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
               const double G,
               const bool compatibleVelocityUpdate);

  //! Destructor.
  virtual ~NBodyGravity();

  //! We augment the generic body force state.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

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

  //! Beginning of timestep work.
  virtual void preStepInitialize(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  //! End of timestep finalizations.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBaseSpace::DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  //! Required label for Physics interface.
  virtual std::string label() const override { return "NBodyGravity"; }

  //! This package opts out of building connectivity.
  virtual bool requireConnectivity() const override { return false; }

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const override;

  //! Return the gravitational potential created by the particle distribution.
  const FieldSpace::FieldList<Dimension, Scalar>& potential() const;

  //! Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The current softening length.
  double softeningLength() const;
  void softeningLength(const double x);

  //! Flag for using the compatible velocity update.
  bool compatibleVelocityUpdate() const;
  void compatibleVelocityUpdate(const bool x);

private:
  
  //! The gravitational potential of the particles.
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential0;
  mutable FieldSpace::FieldList<Dimension, Scalar> mVel02;

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
  
  //! Flag for compatible velocity.
  bool mCompatibleVelocityUpdate;

  // Default constructor -- disabled.
  NBodyGravity();

  // Copy constructor -- disabled.
  NBodyGravity(const NBodyGravity&);

  // Assignment operator -- disabled.
  NBodyGravity& operator=(const NBodyGravity&);

  // Worker for accumulating pair-wise forces.
  void applyPairForces(const std::vector<Scalar>& otherMass,
                       const std::vector<Vector>& otherPosition,
                       const FieldSpace::FieldList<Dimension, Vector>& position,
                       FieldSpace::FieldList<Dimension, Vector>& DvDt,
                       FieldSpace::FieldList<Dimension, Scalar>& potential) const;

  // Methods for serializing/deserializing point values.
  void serialize(const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                 std::vector<char>& buffer) const;
  void deserialize(const std::vector<char>& buffer,
                   std::vector<typename Dimension::Scalar>& mass,
                   std::vector<typename Dimension::Vector>& position) const;

}; // end class NBodyGravity

} // end namespace GravitySpace

} // end namespace Spheral

#endif
