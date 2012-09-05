//---------------------------------Spheral++----------------------------------//
// GadgetGravityForce -- a gravity force that uses the Gadget gravity code.
//
// Created by JNJ, Tue Jul 23 21:06:23 PDT 2002
//----------------------------------------------------------------------------//
#ifndef GadgetGravityForce_HH
#define GadgetGravityForce_HH

#ifdef USE_GADGET

#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

namespace Spheral {
namespace GadgetSpace {

using namespace Spheral::PhysicsSpace;

class GravityForce: public PhysicsSpace::GenericBodyForce<Dim<3> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Constructors.
  GravityForce();

  // Destructor.
  virtual ~GravityForce();

  // This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const DataBase<Dim<3> >& dataBase,
                           const Scalar& time,
                           const Scalar& dt,
                           FieldList<Dim<3>, Vector>& DvDt) const;

  // Vote on Gadget's timestep.
  virtual Scalar dt(const DataBase<Dim<3> >& dataBase, Scalar currentTime) const;

  // Make sure that Gadget's internal state is initialized before cycling.
  virtual void initialize(const DataBase<Dim<3> >& db, 
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                          const Scalar& time, const Scalar& dt);

  // Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  // Return the gravitational potential created by the particle distribution.
  const FieldList<Dim<3> , Scalar>& potential() const;

  // Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

private:
  
  // This static variable keeps track of whether we have initialized a 
  // GadgetGravityForce object.  If we have, then we may not initialize any 
  // others.  Gadget's internal state is given by a bunch of global variables,
  // so we can only have one such object at a time.
  static bool mIsInitialized;

  // The size of the Gadget particle list allocated.
  size_t mSizeOfGadgetParticleList;

  // The gravitational potential of the particles.
  // This is mutable because it is not considered part of the 
  // nominal state of the physics.
  mutable FieldList<Dim<3> , Scalar> mPotential;

  // The total potential energy of the particles.  Also mutable.
  mutable Scalar mExtraEnergy;
  
  // Initialize Gadget.
  void mInitializeGadget(const DataBase<Dim<3> >& db);

  // Copy constructor -- disabled.
  GravityForce(const GravityForce&);

  // Assignment operator -- disabled.
  GravityForce& operator=(const GravityForce&);

}; // end class GravityForce

} // end namespace GadgetSpace

} // end namespace Spheral

#endif // USE_GADGET

#endif
