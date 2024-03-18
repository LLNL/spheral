//---------------------------------Spheral++----------------------------------//
// CRKSPHHydroBase -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHHydroBase_hh__
#define __Spheral_CRKSPHHydroBase_hh__

#include "Physics/GenericHydro.hh"
#include "Geometry/CellFaceFlag.hh"
#include "RK/RKCorrectionParams.hh"

#include <string>

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
}

namespace Spheral {

template<typename Dimension>
class CRKSPHHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  DataBase<Dimension>& dataBase,
                  ArtificialViscosity<Dimension>& Q,
                  const RKOrder order,
                  const double filter,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool XSPH,
                  const MassDensityType densityUpdate,
                  const HEvolutionType HUpdate,
                  const double epsTensile,
                  const double nTensile);

  // Destructor.
  virtual ~CRKSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                            State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivs) override;

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

  // We require RK corrections
  virtual std::set<RKOrder> requireReproducingKernels() const override;

  // The spatial order
  RKOrder correctionOrder() const;
  void correctionOrder(RKOrder val);

  // Flag to choose whether we want to sum for density, or integrate
  // the continuity equation.
  MassDensityType densityUpdate() const;
  void densityUpdate(MassDensityType type);

  // Flag to select how we want to evolve the H tensor.
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

  // The object defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // Fraction of centroidal filtering to apply.
  double filter() const;
  void filter(double val);

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);
    
  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldList<Dimension, Scalar>&    entropy() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldList<Dimension, Scalar>&    effectiveViscousPressure() const;
  const FieldList<Dimension, Scalar>&    viscousWork() const;
  const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldList<Dimension, Vector>&    massFirstMoment() const;
  const FieldList<Dimension, SymTensor>& massSecondMomentEta() const;
  const FieldList<Dimension, SymTensor>& massSecondMomentLab() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;

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
  virtual std::string label() const override { return "CRKSPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  RKOrder mOrder;
  MassDensityType mDensityUpdate;
  HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mEvolveTotalEnergy, mXSPH;
  double mfilter;
  Scalar mEpsTensile, mnTensile;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;
  FieldList<Dimension, Scalar>    mEntropy;

  FieldList<Dimension, SymTensor> mHideal;
  FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldList<Dimension, Scalar>    mEffViscousPressure;
  FieldList<Dimension, Scalar>    mViscousWork;

  FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldList<Dimension, Vector>    mMassFirstMoment;
  FieldList<Dimension, SymTensor> mMassSecondMomentEta;
  FieldList<Dimension, SymTensor> mMassSecondMomentLab;

  FieldList<Dimension, Vector>    mXSPHDeltaV;
  FieldList<Dimension, Vector>    mDxDt;

  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDHDt;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;

  std::vector<Vector>             mPairAccelerations;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  CRKSPHHydroBase();
  CRKSPHHydroBase(const CRKSPHHydroBase&);
  CRKSPHHydroBase& operator=(const CRKSPHHydroBase&);
};

}

#include "CRKSPHHydroBaseInline.hh"

#endif
