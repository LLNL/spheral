//---------------------------------Spheral++----------------------------------//
// FVPMTotalHydroBase -- The FVPM/AFVPM hydrodynamic package for Spheral++.
//
// Created by JMO, Tue Aug  3 21:22:10 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_FVPMTotalHydroBase_hh__
#define __Spheral_FVPMTotalHydroBase_hh__

#include <string>

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
class FVPMTotalHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  FVPMTotalHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                     const TableKernel<Dimension>& W,
                     ArtificialViscosity<Dimension>& Q,
                     const double cfl,
                     const bool useVelocityMagnitudeForDt,
                     const HEvolutionType HUpdate);

  // Destructor.
  virtual ~FVPMTotalHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase);

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs);

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs);
                       
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Finalize the hydro at the completion of an integration step.
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs);
               
  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  HEvolutionType HEvolution() const;
  void HEvolution(const HEvolutionType type);

  // The object defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    positionWeight() const;
  const FieldList<Dimension, Scalar>&    omegaGradh() const;
  const FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldList<Dimension, Scalar>&    massDensitySum() const;
  const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldList<Dimension, std::vector<Vector> >& pairAccelerations() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "FVPMTotalHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  MassDensityType mDensityUpdate;
  HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mGradhCorrection, mXSPH;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mPositionWeight;
  FieldList<Dimension, Scalar>    mOmegaGradh;
  FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldList<Dimension, SymTensor> mHideal;
  FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldList<Dimension, Scalar>    mMassDensitySum;

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

  FieldList<Dimension, std::vector<Vector> > mPairAccelerations;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  FVPMTotalHydroBase();
  FVPMTotalHydroBase(const FVPMTotalHydroBase&);
  FVPMTotalHydroBase& operator=(const FVPMTotalHydroBase&);
};

}

#include "FVPMTotalHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class FVPMTotalHydroBase;
}

#endif
