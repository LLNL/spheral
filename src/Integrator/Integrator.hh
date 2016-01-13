//---------------------------------Spheral++----------------------------------//
// Integrator -- The topmost abstract base class for all integrator classes
// Spheral++.  Integrator classes take a list of Physics packages and advance
// them in time.
//
// Created by JMO, Wed May 31 21:58:08 PDT 2000
//----------------------------------------------------------------------------//
#ifndef Integrator_HH
#define Integrator_HH

#include <string>
#ifndef __GCCXML__
#include <vector>
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace PhysicsSpace {
    template<typename Dimension> class Physics;
  }
  namespace FileIOSpace {
    class FileIO;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
}

namespace Spheral {
namespace IntegratorSpace {

template<typename Dimension>
class Integrator {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<PhysicsSpace::Physics<Dimension>*>::iterator PackageIterator;
  typedef typename std::vector<PhysicsSpace::Physics<Dimension>*>::const_iterator ConstPackageIterator;

  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::iterator BoundaryIterator;
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Constructors.
  Integrator();
  Integrator(DataBaseSpace::DataBase<Dimension>& dataBase);
  Integrator(DataBaseSpace::DataBase<Dimension>& dataBase,
             const std::vector<PhysicsSpace::Physics<Dimension>*>& physicsPackages);

  // Destructor.
  virtual ~Integrator();

  // Assignment.
  Integrator& operator=(const Integrator& rhs);

  // All Integrator classes must define the dt and step methods.
  virtual void step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) = 0;
  virtual void step(Scalar maxTime);

  // Provide a method of looping over the physics packages and picking a
  // time step.
  virtual Scalar selectDt(const Scalar dtMin, 
                          const Scalar dtMax,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs) const;

  // Perform any expected initializations of the integrator.
  virtual void initialize(State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs);

  // Perform generic initializations at the beginning of a timestep.
  virtual void preStepInitialize(const double t,
                                 const double dt,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);

  // Finalize at the end a timestep.
  virtual void finalize(const double t,
                        const double dt,
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Iterate over all physics packages and call evaluateDerivatives.
  void evaluateDerivatives(const Scalar t,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Iterate over all physics packages and call finalizeDerivatives.
  void finalizeDerivatives(const Scalar t,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Iterate over all physics packages and call postStateUpdate
  void postStateUpdate(const DataBaseSpace::DataBase<Dimension>& dataBase,
                       State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivs) const;

  // Add a Physics package.
  void appendPhysicsPackage(PhysicsSpace::Physics<Dimension>& package);

  // Test if the given Physics package is listed in the integrator.
  bool havePhysicsPackage(const PhysicsSpace::Physics<Dimension>& package) const;

  // Get the unique set of boundary conditions across all physics packages.
  std::vector<BoundarySpace::Boundary<Dimension>*> uniqueBoundaryConditions() const;

  // Set the ghost nodes for all node lists according to the boundary 
  // conditions.
  void setGhostNodes();

  // Set the ghost node values on the Fields of the nodes lists in the
  // data base.
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Finalize the ghost node boundary conditions.
  void finalizeGhostBoundaries();

  // Find the nodes in violation of the boundary conditions.
  void setViolationNodes();

  // Rest any internal nodes in violation of boundary conditions to be brought 
  // into compliance.
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  // Copy the ghost positions and H's from one state to another.
  void copyGhostState(const State<Dimension>& state0,
                      State<Dimension>& state1) const;

  // Advance the set of Physics packages to the given time.
  virtual void advance(Scalar goalTime);

  // Access the current time.
  Scalar currentTime() const;
  void currentTime(Scalar time);

  // Access the current cycle.
  int currentCycle() const;
  void currentCycle(int cycle);

  // Access the minimum allowed timestep.
  Scalar dtMin() const;
  void dtMin(Scalar dt);

  // Access the maximum allowed timestep.
  Scalar dtMax() const;
  void dtMax(Scalar dt);

  // Access the last timestep.
  Scalar lastDt() const;
  void lastDt(Scalar dt);

  // Access the timestep growth factor.
  Scalar dtGrowth() const;
  void dtGrowth(Scalar fraction);

  // Public const access to the DataBase.
  const DataBaseSpace::DataBase<Dimension>& dataBase() const;

  // Access the list of physics packages.
  const std::vector<PhysicsSpace::Physics<Dimension>*>& physicsPackages() const;

  // Provide standard iterator methods over the physics package list.
  PackageIterator physicsPackagesBegin();
  PackageIterator physicsPackagesEnd();

  ConstPackageIterator physicsPackagesBegin() const;
  ConstPackageIterator physicsPackagesEnd() const;

  // Flag to determine whether or not to be rigorous about about boundaries.
  bool rigorousBoundaries() const;
  void rigorousBoundaries(bool value);

  // If we're not being rigorous about boundary conditions, how frequently
  // do we update them?
  int updateBoundaryFrequency() const;
  void updateBoundaryFrequency(const int value);

  // Select whether the integrator is verbose or not during a cycle.
  bool verbose() const;
  void verbose(bool x);

  // Timestep threshold where we start reporting the timestep reasons even if not verbose.
  Scalar dtThreshold() const;
  void dtThreshold(const Scalar x);

  // Select whether we should run in a mode the ensures domain decomposition independence.
  // Possibly some performance impact.
  bool domainDecompositionIndependent() const;
  void domainDecompositionIndependent(const bool x);

  // Select whether we're going to enforce culling of ghost nodes or not.
  bool cullGhostNodes() const;
  void cullGhostNodes(const bool x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "Integrator"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //-------------------------- Protected Interface --------------------------//
  // Allow write access to the DataBase for descendent classes.
  DataBaseSpace::DataBase<Dimension>& accessDataBase();

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mDtMin, mDtMax, mDtGrowth, mLastDt, mCurrentTime, mDtThreshold;
  int mCurrentCycle, mUpdateBoundaryFrequency;
  bool mVerbose, mRequireConnectivity, mRequireGhostConnectivity;
  DataBaseSpace::DataBase<Dimension>* mDataBasePtr;
  std::vector<PhysicsSpace::Physics<Dimension>*> mPhysicsPackages;
  bool mRigorousBoundaries, mCullGhostNodes;

#ifndef __GCCXML__
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif
};

}
}

#ifndef __GCCXML__
#include "IntegratorInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace IntegratorSpace {
    template<typename Dimension> class Integrator;
  }
}

#endif
