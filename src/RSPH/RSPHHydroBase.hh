//---------------------------------Spheral++----------------------------------//
// RSPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_RSPHHydroBase_hh__
#define __Spheral_RSPHHydroBase_hh__

#include <string>

#include "Physics/Physics.hh"
#include "Physics/GenericHydro.hh"
namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class SmoothingScaleBase;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class RSPHHydroBase: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  RSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
               DataBase<Dimension>& dataBase,
               const TableKernel<Dimension>& W,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool evolveTotalEnergy,
               const bool XSPH,
               const bool correctVelocityGradient,
               const HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax);

  // Destructor.
  virtual ~RSPHHydroBase();

    // We require all Physics packages to provide a method returning their vote
  // for the next time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  
  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                       
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;


  //Evaluate Derivatives sub--routines
  void
  evaluateSpatialGradients(const typename Dimension::Scalar time,
                         const typename Dimension::Scalar dt,
                         const DataBase<Dimension>& dataBase,
                         const State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives) const;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;


  void computeHLLCstate( const Vector& rij,
                         int nodeListi,
                         int nodeListj,
                         int i,
                         int j,
                         const Scalar& Pi,  
                         const Scalar& Pj,
                         const Scalar& rhoi, 
                         const Scalar& rhoj,
                         const Vector& vi,   
                         const Vector& vj,
                         const Scalar& ci,   
                         const Scalar& cj, 
                         const Vector& DpDxi,
                         const Vector& DpDxj,
                         const Tensor& DvDxi,
                         const Tensor& DvDxj,
                         Vector& vstar,
                         Scalar& Pstar) const;

  // Also allow access to the CFL timestep safety criteria.
  Scalar cfl() const;
  void cfl(Scalar cfl);

  // Attribute to determine if the absolute magnitude of the velocity should
  // be used in determining the timestep.
  bool useVelocityMagnitudeForDt() const;
  void useVelocityMagnitudeForDt(bool x);

  // Return the cumulative neighboring statistics.
  int minMasterNeighbor() const;
  int maxMasterNeighbor() const;
  double averageMasterNeighbor() const;

  int minCoarseNeighbor() const;
  int maxCoarseNeighbor() const;
  double averageCoarseNeighbor() const;

  int minRefineNeighbor() const;
  int maxRefineNeighbor() const;
  double averageRefineNeighbor() const;

  int minActualNeighbor() const;
  int maxActualNeighbor() const;
  double averageActualNeighbor() const;

  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  HEvolutionType HEvolution() const;
  void HEvolution(HEvolutionType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(bool val);

  // Flag controlling if we evolve total or specific energy.
  bool evolveTotalEnergy() const;
  void evolveTotalEnergy(bool val);

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH() const;
  void XSPH(bool val);

  // Flag to determine if we're applying the linear correction for the velocity gradient.
  bool correctVelocityGradient() const;
  void correctVelocityGradient(bool val);

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // Access the stored interpolation kernels.
  const TableKernel<Dimension>& kernel() const;

  // The object defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, Scalar>&    normalization() const;
  const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Tensor>&    M() const;
  const FieldList<Dimension, Tensor>&    localM() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const std::vector<Vector>&             pairAccelerations() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "RSPHHydroBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //-------------------------- Protected Interface --------------------------//
  void updateMasterNeighborStats(int numMaster) const;
  void updateCoarseNeighborStats(int numNeighbor) const;
  void updateRefineNeighborStats(int numNeighbor) const;
  void updateActualNeighborStats(int numNeighbor) const;


  // The interpolation kernels.
  const TableKernel<Dimension>& mKernel;

  // The method defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mEvolveTotalEnergy,  mXSPH, mCorrectVelocityGradient;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;

  FieldList<Dimension, SymTensor> mHideal;
  FieldList<Dimension, Scalar>    mNormalization;

  FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldList<Dimension, Vector>    mXSPHDeltaV;

  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDHDt;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;
  FieldList<Dimension, Tensor>    mM;
  FieldList<Dimension, Tensor>    mLocalM;

  std::vector<Vector>             mPairAccelerations;

  FieldList<Dimension, Vector> mDpDx;
  FieldList<Dimension, Tensor> mLastDvDx;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;

private:

  Scalar mCfl;
  bool mUseVelocityMagnitudeForDt;

  mutable int mMinMasterNeighbor, mMaxMasterNeighbor, mSumMasterNeighbor;
  mutable int mMinCoarseNeighbor, mMaxCoarseNeighbor, mSumCoarseNeighbor;
  mutable int mMinRefineNeighbor, mMaxRefineNeighbor, mSumRefineNeighbor;
  mutable int mMinActualNeighbor, mMaxActualNeighbor, mSumActualNeighbor;
  mutable int mNormMasterNeighbor;
  mutable int mNormCoarseNeighbor;
  mutable int mNormRefineNeighbor;
  mutable int mNormActualNeighbor;
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  RSPHHydroBase();
  RSPHHydroBase(const RSPHHydroBase&);
  RSPHHydroBase& operator=(const RSPHHydroBase&);
};

}

#include "RSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RSPHHydroBase;
}

#endif
