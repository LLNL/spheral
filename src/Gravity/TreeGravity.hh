//---------------------------------Spheral++----------------------------------//
// TreeGravity -- An implementation of the Barnes-Hut tree n-body gravity
//                solver.
//
//! \author $Author: mikeowen $
//! \version $Revision: 3174 $
//! \date $Date: 2009-07-22 14:23:24 -0700 (Wed, 22 Jul 2009) $
//
//----------------------------------------------------------------------------//
#ifndef __Spheral_TreeGravity__
#define __Spheral_TreeGravity__

#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace GravitySpace {

template <typename Dimension>
class TreeGravity: public PhysicsSpace::GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;

  //! Constructor.
  //! \param G -- the gravitational constant.
  //! \param opening -- the opening ratio for approximating forces.
  TreeGravity(const double G,
              const double opening);

  //! Destructor.
  virtual ~TreeGravity();

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

  //! The opening angle threshold when we shift to tree cell approximations.
  double opening() const;

private:
  double mG, mOpening;

  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;
  mutable Scalar mExtraEnergy;
  
  // Default constructor -- disabled.
  TreeGravity();

  // Copy constructor -- disabled.
  TreeGravity(const TreeGravity&);

  // Assignment operator -- disabled.
  TreeGravity& operator=(const TreeGravity&);

};

}
}

#endif
