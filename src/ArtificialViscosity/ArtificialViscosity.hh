//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//
// Created by JMO, Sun May 21 21:16:43 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosity__
#define __Spheral_ArtificialViscosity__

#include "Physics/Physics.hh"
#include "Field/FieldList.hh"
#include "DataOutput/registerWithRestart.hh"

#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class Boundary;
class FileIO;

template<typename Dimension, typename QPiType>
class ArtificialViscosity: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ReturnType = QPiType;

  // Constructors, destructor
  ArtificialViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& kernel);
  virtual ~ArtificialViscosity() = default;

  // No default constructor, copying, or assignment
  ArtificialViscosity() = delete;
  ArtificialViscosity(const ArtificialViscosity&) = delete;
  ArtificialViscosity& operator=(const ArtificialViscosity&) = delete;

  //...........................................................................
  // Virtual methods we expect ArtificialViscosities to provide
  // Some AVs need the velocity gradient computed, so they should override this to true
  virtual bool requireVelocityGradient()                                  const { return false; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
  virtual void QPiij(QPiType& QPiij, QPiType& QPiji,    // result for QPi (Q/rho^2)
                     Scalar& Qij, Scalar& Qji,          // result for viscous pressure
                     const unsigned nodeListi, const unsigned i, 
                     const unsigned nodeListj, const unsigned j,
                     const Vector& xi,
                     const SymTensor& Hi,
                     const Vector& etai,
                     const Vector& vi,
                     const Scalar rhoi,
                     const Scalar csi,
                     const Vector& xj,
                     const SymTensor& Hj,
                     const Vector& etaj,
                     const Vector& vj,
                     const Scalar rhoj,
                     const Scalar csj,
                     const FieldList<Dimension, Scalar>& fCl,
                     const FieldList<Dimension, Scalar>& fCq,
                     const FieldList<Dimension, Tensor>& DvDx) const = 0;

  //...........................................................................
  // Standard Physics package methods
  // Most ArtificialViscosities will not have an evaluateDerivatives, so by default no-op this
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override {};

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  // Post-state update is our chance to update the velocity gradient if needed
  virtual 
  bool postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  // Update the locally stored velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs);

  // Calculate the curl of the velocity given the stress tensor.
  Scalar curlVelocityMagnitude(const Tensor& DvDx) const;

  // Find the Balsara shear correction multiplier
  Scalar calcBalsaraShearCorrection(const Tensor& DvDx,
                                    const SymTensor& H,
                                    const Scalar& cs) const;

  // Access stored state
  Scalar                              Cl()                                const { return mClinear; }
  Scalar                              Cq()                                const { return mCquadratic; }
  bool                                balsaraShearCorrection()            const { return mBalsaraShearCorrection; }
  Scalar epsilon2()                                                       const { return mEpsilon2; }
  Scalar negligibleSoundSpeed()                                           const { return mNegligibleSoundSpeed; }
  const FieldList<Dimension, Scalar>& maxViscousPressure()                const { return mMaxViscousPressure; }
  const FieldList<Dimension, Tensor>& DvDx()                              const { return mDvDx; }

  bool                                rigorousVelocityGradient()          const { return mRigorousVelocityGradient; }
  const TableKernel<Dimension>&       kernel()                            const { return mWT; }

  void Cl(Scalar x)                                                             { mClinear = x; }
  void Cq(Scalar x)                                                             { mCquadratic = x; }
  void balsaraShearCorrection(bool x)                                           { mBalsaraShearCorrection = x; }
  void epsilon2(Scalar x)                                                       { mEpsilon2 = x; }
  void negligibleSoundSpeed(Scalar x)                                           { REQUIRE(x > 0.0); mNegligibleSoundSpeed = x; }
  void rigorousVelocityGradient(bool x)                                         { mRigorousVelocityGradient = x; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "ArtificialViscosity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  Scalar mClinear;
  Scalar mCquadratic;

  // Switch for the Balsara shear correction.
  bool mBalsaraShearCorrection;

  // Parameters for the Q limiter.
  Scalar mEpsilon2;
  Scalar mNegligibleSoundSpeed;
    
  // Maintain the last max viscous pressure for timestep control
  FieldList<Dimension, Scalar> mMaxViscousPressure;

  // State for maintaining the velocity gradient
  bool mRigorousVelocityGradient;
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Tensor> mM;
  FieldList<Dimension, Tensor> mDvDx;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "ArtificialViscosityInline.hh"

#endif
