//---------------------------------Spheral++----------------------------------//
// Physics -- The topmost abstract base class for all physics packages in
// Spheral++.  Physics defines the methods generic to all physics classes,
// but leaves specialized interfaces for particular types of physics (radiative vs.
// hydro, etc.) to descendent interface classes.
//
// Created by JMO, Wed May 24 10:26:36 PDT 2000
//----------------------------------------------------------------------------//
#ifndef Physics_HH
#define Physics_HH

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include <string>

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
}

namespace Spheral {
namespace PhysicsSpace {

template<typename Dimension>
class Physics {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::iterator BoundaryIterator;
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;
  typedef typename std::pair<double, std::string> TimeStepType;

  // Constructors.
  Physics();

  // Destructor.
  virtual ~Physics();

  //******************************************************************************//
  // Methods all Physics packages must provide.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const = 0;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const = 0;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state) = 0;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) = 0;

  //******************************************************************************//
  // Methods for handling boundary conditions.
  // Add a Boundary condition.
  void appendBoundary(BoundarySpace::Boundary<Dimension>& boundary);  // To end of boundary list
  void prependBoundary(BoundarySpace::Boundary<Dimension>& boundary); // To beginning of boundary list
  void clearBoundaries();                                             // Remove all boundary conditions

  // Test if the given Boundary condition is registered.
  bool haveBoundary(const BoundarySpace::Boundary<Dimension>& boundary) const;

  // Provide standard iterator methods over the boundary conditions list.
  BoundaryIterator boundaryBegin();
  BoundaryIterator boundaryEnd();

  ConstBoundaryIterator boundaryBegin() const;
  ConstBoundaryIterator boundaryEnd() const;

  // Access the list of boundary conditions.
  const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaryConditions() const;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);

  //******************************************************************************//
  // An optional hook to initialize once when the problem is starting up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Some packages might want a hook to do some initializations before the
  // evaluateDerivatives() method is called.
  virtual void initialize(const Scalar time, 
                          const Scalar dt,
                          const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs);

  // Similarly packages might want a hook to do some post-step finalizations.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBaseSpace::DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Provide a hook to be called after all physics packages have had their
  // evaluateDerivatives method called, but before anyone does anything
  // with those derivatives.
  virtual 
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Provide a hook to be called after the state has been updated and 
  // boundary conditions have been enforced.
  virtual 
  void postStateUpdate(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivatives) const;

  // Some physics does not require the connectivity be constructed.
  virtual bool requireConnectivity() const;

  // Many physics packages will have their own representations of energy in the
  // system (gravitational potential energy, radiative losses, etc.)
  virtual Scalar extraEnergy() const;

  // Many physics packages will also have their own representations of momentum in the
  // system (electromagnetic momentum flux density, etc.) 
  virtual Vector extraMomentum() const;

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  std::vector<BoundarySpace::Boundary<Dimension>*> mBoundaryConditions;
#endif
};

}
}

#ifndef __GCCXML__
#include "PhysicsInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class Physics;
  }
}

#endif
