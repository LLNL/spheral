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
namespace FVPMSpace {

template<typename Dimension>
class FVPMTotalHydroBase: public PhysicsSpace::GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  FVPMTotalHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                     const KernelSpace::TableKernel<Dimension>& W,
                     ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                     const double cfl,
                     const bool useVelocityMagnitudeForDt,
                     const HEvolutionType HUpdate);

  // Destructor.
  virtual ~FVPMTotalHydroBase();

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
  // the continuity equation.
  HEvolutionType HEvolution() const;
  void HEvolution(const HEvolutionType type);

  // The object defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, int>&       timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    positionWeight() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    omegaGradh() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    massDensitySum() const;
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

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "FVPMTotalHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  MassDensityType mDensityUpdate;
  HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mGradhCorrection, mXSPH;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, int>       mTimeStepMask;
  FieldSpace::FieldList<Dimension, Scalar>    mPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldSpace::FieldList<Dimension, Scalar>    mPositionWeight;
  FieldSpace::FieldList<Dimension, Scalar>    mOmegaGradh;
  FieldSpace::FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  FieldSpace::FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mMassDensitySum;

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

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  FVPMTotalHydroBase();
  FVPMTotalHydroBase(const FVPMTotalHydroBase&);
  FVPMTotalHydroBase& operator=(const FVPMTotalHydroBase&);
};

}
}

#ifndef __GCCXML__
#include "FVPMTotalHydroBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace FVPMSpace {
    template<typename Dimension> class FVPMTotalHydroBase;
  }
}

#endif
