//---------------------------------Spheral++----------------------------------//
// TaylorSPHHydroBase -- The TaylorSPH hydrodynamic package for Spheral++.
// Based on G. Oger et al., JCP 225, pp 1472-1492 (2007)
//
// Created by JMO, Mon Jun 30 21:36:40 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_TaylorSPHHydroBase_hh__
#define __Spheral_TaylorSPHHydroBase_hh__

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
class TaylorSPHHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  TaylorSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                const TableKernel<Dimension>& W,
                ArtificialViscosity<Dimension>& Q,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool XSPH,
                const HEvolutionType HUpdate);

  // Destructor.
  virtual ~TaylorSPHHydroBase();

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
  HEvolutionType HEvolution() const;
  void HEvolution(const HEvolutionType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(const bool val);

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH() const;
  void XSPH(const bool val);

  // The object defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
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

  const FieldList<Dimension, Tensor>&    D() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "TaylorSPHHydroBase"; }
  virtual void dumpState(FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIO& file, std::string pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mXSPH;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldList<Dimension, SymTensor> mHideal;
  FieldList<Dimension, Scalar>    mMaxViscousPressure;

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

  FieldList<Dimension, Tensor>    mD;

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  TaylorSPHHydroBase();
  TaylorSPHHydroBase(const TaylorSPHHydroBase&);
  TaylorSPHHydroBase& operator=(const TaylorSPHHydroBase&);
};

}

#include "TaylorSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class TaylorSPHHydroBase;
}

#endif
