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

#include "RK/RKCorrectionParams.hh"

#include <vector>
#include <set>
#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;

template<typename Dimension>
class Physics {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;
  using TimeStepType = typename std::pair<double, std::string>;

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
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const = 0;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const = 0;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) = 0;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) = 0;

  // It's useful to have labels for Physics packages.  We'll require this to have
  // the same signature as the restart label.
  virtual std::string label() const = 0;

  //******************************************************************************//
  // Methods for handling boundary conditions.
  // Add a Boundary condition.
  void appendBoundary(Boundary<Dimension>& boundary);  // To end of boundary list
  void prependBoundary(Boundary<Dimension>& boundary); // To beginning of boundary list
  void clearBoundaries();                              // Remove all boundary conditions

  // Test if the given Boundary condition is registered.
  bool haveBoundary(const Boundary<Dimension>& boundary) const;

  // Provide standard iterator methods over the boundary conditions list.
  BoundaryIterator boundaryBegin();
  BoundaryIterator boundaryEnd();

  ConstBoundaryIterator boundaryBegin() const;
  ConstBoundaryIterator boundaryEnd() const;

  // Access the list of boundary conditions.
  const std::vector<Boundary<Dimension>*>& boundaryConditions() const;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);

  //******************************************************************************//
  // The following two methods are for adding additional physics package that
  // should be run before or after this one.  These methods are intended for use
  // in our Python interface where we can automatically construct packages for the
  // user for simplicities sake.
  // Add optional Physics packages that should be inserted in the Physics package list after this one
  void appendSubPackage(Physics<Dimension>& package);

  // Add optional Physics packages that should be inserted in the Physics package list before this one
  void prependSubPackage(Physics<Dimension>& package);

  // Access the sets of pre- and post-subpackages
  const std::vector<Physics<Dimension>*>& postSubPackages() const;
  const std::vector<Physics<Dimension>*>& preSubPackages() const;

  //******************************************************************************//
  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase);

  // A second optional method to be called on startup, after Physics::initializeProblemStartup
  // has been called.
  // This method is called after independent variables have been initialized and put into
  // the state and derivatives. During this method, the dependent state, such as
  // temperature and pressure, is initialized so that all the fields in the initial
  // state and derivatives objects are valid.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs);

  // Optional hook to be called at the beginning of a time step.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);

  // Some packages might want a hook to do some initializations before the
  // evaluateDerivatives() method is called.
  virtual void initialize(const Scalar time, 
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs);

  // Similarly packages might want a hook to do some post-step finalizations.
  // Really we should rename this post-step finalize.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Provide a hook to be called after all physics packages have had their
  // evaluateDerivatives method called, but before anyone does anything
  // with those derivatives.
  virtual 
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Provide a hook to be called after the state has been updated and 
  // boundary conditions have been enforced.
  virtual 
  void postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives);

  // Some physics does not require the connectivity be constructed.
  virtual bool requireConnectivity() const;

  // Some physics algorithms require ghost connectivity to be constructed.
  virtual bool requireGhostConnectivity() const;

  // Some physics algorithms require overlap connectivity.
  virtual bool requireOverlapConnectivity() const;

  // Some physics algorithms require intersection connectivity
  virtual bool requireIntersectionConnectivity() const;

  // Does this package require reproducing kernel functions?
  virtual std::set<RKOrder> requireReproducingKernels() const;

  // If using reproducing kernels, do we need the second derivative?
  virtual bool requireReproducingKernelHessian() const;

  // Does this package need an update of reproducing kernels during finalize?
  virtual bool updateReproducingKernelsInFinalize() const;
  
  // Many physics packages will have their own representations of energy in the
  // system (gravitational potential energy, radiative losses, etc.)
  virtual Scalar extraEnergy() const;

  // Many physics packages will also have their own representations of momentum in the
  // system (electromagnetic momentum flux density, etc.) 
  virtual Vector extraMomentum() const;

  // Register any additional state for visualization.
  virtual void registerAdditionalVisualizationState(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state);

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<Boundary<Dimension>*> mBoundaryConditions;
  std::vector<Physics<Dimension>*> mPreSubPackages, mPostSubPackages;
};

}

#include "PhysicsInline.hh"

#endif
