//---------------------------------Spheral++----------------------------------//
// MHD -- Resistive MHD done SPH style.  This object cooperates with the 
//        SPH Hydro class.
//
//! \author $Author: jeffjohnson $
//! \version $Revision: 2239 $
//! \date $Date: 2007-05-28 23:58:39 -0700 (Mon, 28 May 2007) $
//
//----------------------------------------------------------------------------//
#ifndef MHD_HH
#define MHD_HH

#include "Python.h"
#include "Physics/Physics.hh"
#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "boost/python/handle.hpp"
#ifndef __GCCXML__
#include "petsc.h"
#include "petscksp.h"
#include <map>
#endif

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

class MHD: public Physics<Dim<3> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  typedef Physics<Dim<3> >::TimeStepType TimeStepType;
  typedef Physics<Dim<3> >::ConstBoundaryIterator ConstBoundaryIterator;

  // Options for evolving the specific thermal energy.
  enum SpecificThermalEnergyType {
    dontUpdateSpecificThermalEnergy = 0,
    integrateSpecificThermalEnergy = 1,
    computeCompatibleSpecificThermalEnergy = 2
  };

  // Options for evolving the total specific energy.
  enum TotalSpecificEnergyType {
    dontUpdateTotalSpecificEnergy = 0,
    integrateTotalSpecificEnergy = 1,
  };

  // Options for dealing with the MHD tensile instability.
  enum TensileStabilizerType {
    noStabilizer = 0,
    maxStressStabilizer = 1,
    MorrisStabilizer = 2,
    B0rveStabilizer = 3,
    extFieldStabilizer = 4
  };
  
  // Options for dealing with the divergence of B.
  enum BDivergenceCleanerType {
    noCleaner = 0,
    GreensFnProjCleaner = 1,
    BiotSavartProjCleaner = 2,
    hyperbolicCleaner = 3
  };

  //! Constructor.
  //! \param W The kernel used to compute SPH interpolants.
  //! \param mu0 The permeability of free space in vacuum.
  //! \param implicitness The degree of implicitness to use (0 to 1.0) in treating magnetic diffusion.  Defaults to 0.5.
  MHD(const TableKernel<Dim<3> >& W,
      double mu0,
      double implicitness = 0.5);

  //! Destructor.
  virtual ~MHD();

  //! This is the derivative method that all Physics classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dim<3> >& dataBase,
                           const State<Dim<3> >& state,
                           StateDerivatives<Dim<3> >& derivs) const;

  //! Register the state for MHD.
  virtual 
  void registerState(DataBase<Dim<3> >& dataBase,
                     State<Dim<3> >& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dim<3> >& dataBase,
                           StateDerivatives<Dim<3> >& derivs);

  // Apply boundary conditions to the fields.
  virtual void applyGhostBoundaries(State<Dim<3> >& state,
                                    StateDerivatives<Dim<3> >& derivs);

  // Enforce boundary conditions for the fields.
  virtual void enforceBoundaries(State<Dim<3> >& state,
                                 StateDerivatives<Dim<3> >& derivs);

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBase<Dim<3> >& dataBase, 
                          const State<Dim<3> >& state,
                          const StateDerivatives<Dim<3> >& derivs,
                          const Scalar currentTime) const;

  //! Pre-step work.
  virtual void initialize(const Scalar& time, 
                          const Scalar& dt,
                          const DataBase<Dim<3> >& db, 
                          State<Dim<3> >& state,
                          StateDerivatives<Dim<3> >& derivs);

  //! Post-step work.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dim<3> >& db, 
                        State<Dim<3> >& state,
                        StateDerivatives<Dim<3> >& derivs);

  // This method performs work after each stage of integration.  In particular, 
  // any magnetic divergence cleaning we do goes here. 
  void postStateUpdate(const DataBase<Dim<3> >& db, 
                       State<Dim<3> >& state,
                       const StateDerivatives<Dim<3> >& derivatives) const;

  //! Return the magnetic diffusion matrix for inspection.
  boost::python::handle<PyObject> diffusionMatrix() const;

  //! Return the right-hand-side vector of the magnetic diffusion
  //! equation for inspection.
  boost::python::handle<PyObject> diffusionRHS() const;

  //! Return the norm of the residual vector after the most recently converged
  //! linear solve of the magnetic diffusion equation.
  double diffusionResidualNorm() const { return mDiffResNorm; }

  //! Return the number of iterations taken in the most recent linear solve
  //! of the magnetic diffusion equation.
  int numDiffusionIterations() const { return mDiffNumIters; }

  //! Return the permeability of free space in vacuum.
  double mu0() const { return mMu0; }

  //! Return the implicitness of the magnetic diffusion solver.
  double implicitness() const { return mImplicitness; }
  void implicitness(double alpha);

  //! This returns the magnetic field energy.
  virtual Scalar extraEnergy() const;

  //! Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

  //! An enum representing the stabilization method.
  TensileStabilizerType stabilizer() const { return mStabilizer; }
  void stabilizer(TensileStabilizerType stab) { mStabilizer = stab; }

  //! The maximum Maxwell stress found in the node distribution.
  const SymTensor& S0() const { return mS0; }

  //! Any external magnetic field to be used with the external field 
  //! stabilizer.
  const Vector& Bext() const { return mBext; }
  void Bext(const Vector& B) { mBext = B; }

  //! An enum representing the method of updating the thermal energy.
  SpecificThermalEnergyType specificThermalEnergyUpdate() const { return mSpecificThermalEnergyUpdate; }
  void specificThermalEnergyUpdate(SpecificThermalEnergyType update) { mSpecificThermalEnergyUpdate = update; }

  //! An enum representing the method of updating the total specific energy.
  TotalSpecificEnergyType totalSpecificEnergyUpdate() const { return mTotalSpecificEnergyUpdate; }
  void totalSpecificEnergyUpdate(TotalSpecificEnergyType update) { mTotalSpecificEnergyUpdate = update; }

  //! An enum representing the magnetic divergence cleaning method.
  BDivergenceCleanerType divBCleaner() const { return mDivBCleaner; }
  void divBCleaner(BDivergenceCleanerType cleaner) { mDivBCleaner = cleaner; }

  //! The maximum, minimum, and averaged measured divergences of B.
  double maxDivB() const { return mMaxDivB; }
  double minDivB() const { return mMinDivB; }
  double averageDivB() const { return mAvgDivB; }

private:
  
  //! Information about the most recent linear solve.  This stuff is outside 
  //! the GCCXML block because the associated accessors are inlined (so 
  //!  GCCXML has to see it).
  mutable double mDiffResNorm;
  mutable int mDiffNumIters;

  //! Permeability of free space in vacuum.
  double mMu0;

  //! Implicitness parameter (0 <= implicitness <= 1).
  double mImplicitness;

#ifndef __GCCXML__
  //! SPH kernel.
  const TableKernel<Dim<3> >& mKernel;

  //! Vector factory.
  mutable PyObject* mVecFactory;

  //! Matrix factory.
  mutable PyObject* mMatFactory;

  //! Diffusion matrix. 
  mutable PyObject* mDiffMatrix;

  //! Diffusion Right-hand-side vector.
  mutable PyObject* mDiffRHS;

  //! Map of boundary nodes that are related to periodic or distributed 
  //! boundary conditions.
  mutable std::map<std::pair<int, int>, int> mIsPeriodicOrDistributedNode;

  //! PETSc linear solver for magnetic diffusion.
  mutable KSP mDiffSolver;

  // A global indexing scheme for all the nodes in our problem of interest.
  mutable FieldList<Dim<3>, int> mNodeIndices;

  // A mapping of nodes to other nodes that overlap them.
  mutable std::map<int, int> mOverlapNodes;

  // The magnetic field energy.
  mutable double mMagneticEnergy;

  // Compute the diffusion matrix structure given a set of nodal data.
  void mComputeMatrixStructure(const DataBase<Dim<3> >& dataBase,
                               const State<Dim<3> >& state) const;
 
  // Compute the curl-curl operator matrix.
  Mat mCurlCurlMatrix(const DataBase<Dim<3> >& dataBase,
                      const State<Dim<3> >& state) const;
 
  // Apply boundary conditions to the linear system representing the 
  // magnetic diffusion equation.
  void mApplyDiffusionBCs(const DataBase<Dim<3> >& dataBase,
                          const State<Dim<3> >& state, 
                          PyObject* matrix, 
                          PyObject* RHS) const;
  
  // Compute the time derivatives resulting from magnetic induction 
  // add them to the given set of derivatives.
  void mAddDiffusionDerivatives(const DataBase<Dim<3> >& dataBase,
                                const State<Dim<3> >& state,
                                const Scalar time,
                                const Scalar dt,
                                StateDerivatives<Dim<3> >& derivs) const;

  // Compute the divergence of B.
  void mComputeDivB(const DataBase<Dim<3> >& dataBase,
                    State<Dim<3> >& state) const;

  // Clean the divergence of B.
  void mCleanDivB(const DataBase<Dim<3> >& dataBase,
                  State<Dim<3> >& state) const;

#endif

  // The method by which we stabilize the MHD tensile instability.
  TensileStabilizerType mStabilizer;

  // The maximum Maxwell stress found within the nodes.
  mutable SymTensor mS0;

  // The external magnetic field (if any) that is subtracted from the 
  // Maxwell stress tensor to provide "tensile" stability.
  Vector mBext;

  // The divergence cleaner.
  BDivergenceCleanerType mDivBCleaner;

  // The methods of updating the thermal and total specific energy.
  SpecificThermalEnergyType mSpecificThermalEnergyUpdate;
  TotalSpecificEnergyType mTotalSpecificEnergyUpdate;

  // The potential involved in the hyperbolic divergence cleaning scheme and 
  // its time derivative.
  FieldList<Dim<3>, Scalar> mPsi, mDpsiDt;

  // Maximum, minimum, and average measured divergence of B.
  mutable double mMaxDivB, mMinDivB, mAvgDivB;

  mutable bool mFirstStep;

  // Default constructor -- disabled.
  MHD();

  // Copy constructor -- disabled.
  MHD(const MHD&);

  // Assignment operator -- disabled.
  MHD& operator=(const MHD&);

}; // end class MHD

} // end namespace Spheral

#endif
