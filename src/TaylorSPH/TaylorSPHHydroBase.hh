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
  namespace NodeSpace {
    template<typename Dimension> class SmoothingScaleBase;
  }
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace TaylorSPHSpace {

template<typename Dimension>
class TaylorSPHHydroBase: public PhysicsSpace::GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  TaylorSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                const KernelSpace::TableKernel<Dimension>& W,
                ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool XSPH,
                const PhysicsSpace::HEvolutionType HUpdate);

  // Destructor.
  virtual ~TaylorSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs);

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBaseSpace::DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs);
                          
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Finalize the hydro at the completion of an integration step.
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBaseSpace::DataBase<Dimension>& dataBase,
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
  PhysicsSpace::HEvolutionType HEvolution() const;
  void HEvolution(const PhysicsSpace::HEvolutionType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(const bool val);

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH() const;
  void XSPH(const bool val);

  // The object defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, int>&       timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldSpace::FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DxDt() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DvDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldSpace::FieldList<Dimension, std::vector<Vector> >& pairAccelerations() const;

  const FieldSpace::FieldList<Dimension, Tensor>&    D() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "TaylorSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, std::string pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  PhysicsSpace::HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mXSPH;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, int>       mTimeStepMask;
  FieldSpace::FieldList<Dimension, Scalar>    mPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldSpace::FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  FieldSpace::FieldList<Dimension, Scalar>    mMaxViscousPressure;

  FieldSpace::FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldSpace::FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldSpace::FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldSpace::FieldList<Dimension, Vector>    mXSPHDeltaV;

  FieldSpace::FieldList<Dimension, Vector>    mDxDt;
  FieldSpace::FieldList<Dimension, Vector>    mDvDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldSpace::FieldList<Dimension, SymTensor> mDHDt;
  FieldSpace::FieldList<Dimension, Tensor>    mDvDx;
  FieldSpace::FieldList<Dimension, Tensor>    mInternalDvDx;

  FieldSpace::FieldList<Dimension, std::vector<Vector> > mPairAccelerations;

  FieldSpace::FieldList<Dimension, Tensor>    mD;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  TaylorSPHHydroBase();
  TaylorSPHHydroBase(const TaylorSPHHydroBase&);
  TaylorSPHHydroBase& operator=(const TaylorSPHHydroBase&);
};

}
}

#include "TaylorSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace TaylorSPHSpace {
    template<typename Dimension> class TaylorSPHHydroBase;
  }
}

#endif
