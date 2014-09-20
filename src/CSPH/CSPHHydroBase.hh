//---------------------------------Spheral++----------------------------------//
// CSPHHydroBase -- The CSPH/ACSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_CSPHHydroBase_hh__
#define __Spheral_CSPHHydroBase_hh__

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
namespace CSPHSpace {

template<typename Dimension>
class CSPHHydroBase: public PhysicsSpace::GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                const KernelSpace::TableKernel<Dimension>& W,
                const KernelSpace::TableKernel<Dimension>& WPi,
                ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool XSPH,
                const PhysicsSpace::MassDensityType densityUpdate,
                const PhysicsSpace::HEvolutionType HUpdate);

  // Destructor.
  virtual ~CSPHHydroBase();

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

  // Hook called after the state has been updated and boundary conditions have been enforced.
  virtual 
  void postStateUpdate(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivatives) const;

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

  // // We need ghost connectivity to be computed.
  // virtual bool requireGhostConnectivity() const { return true; }

  // Flag to choose whether we want to sum for density, or integrate
  // the continuity equation.
  PhysicsSpace::MassDensityType densityUpdate() const;
  void densityUpdate(const PhysicsSpace::MassDensityType type);

  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  PhysicsSpace::HEvolutionType HEvolution() const;
  void HEvolution(const PhysicsSpace::HEvolutionType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(const bool val);

  // Flag to determine if we're using the grad h correction.
  bool gradhCorrection() const;
  void gradhCorrection(const bool val);

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH() const;
  void XSPH(const bool val);

  // The object defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // Fraction of centroidal filtering to apply.
  double filter() const;
  void filter(const double val);

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, int>&       timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    volume() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  // const FieldSpace::FieldList<Dimension, Scalar>&    massDensitySum() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldSpace::FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DxDt() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DvDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DmassDensityDx() const;
  const FieldSpace::FieldList<Dimension, std::vector<Vector> >& pairAccelerations() const;

  const FieldSpace::FieldList<Dimension, Scalar>&    m0() const;
  const FieldSpace::FieldList<Dimension, Vector>&    m1() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& m2() const;

  const FieldSpace::FieldList<Dimension, Scalar>&    A0() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    A() const;
  const FieldSpace::FieldList<Dimension, Vector>&    B() const;
  const FieldSpace::FieldList<Dimension, Vector>&    C() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    D() const;
  const FieldSpace::FieldList<Dimension, Vector>&    gradA0() const;
  const FieldSpace::FieldList<Dimension, Vector>&    gradA() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    gradB() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "CSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, std::string pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  PhysicsSpace::MassDensityType mDensityUpdate;
  PhysicsSpace::HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mGradhCorrection, mXSPH;
  double mfilter;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, int>       mTimeStepMask;
  FieldSpace::FieldList<Dimension, Scalar>    mPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldSpace::FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  FieldSpace::FieldList<Dimension, Scalar>    mMaxViscousPressure;
  // FieldSpace::FieldList<Dimension, Scalar>    mMassDensitySum;

  FieldSpace::FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldSpace::FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldSpace::FieldList<Dimension, Vector>    mXSPHDeltaV;

  FieldSpace::FieldList<Dimension, Vector>    mDxDt;
  FieldSpace::FieldList<Dimension, Vector>    mDvDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldSpace::FieldList<Dimension, SymTensor> mDHDt;
  FieldSpace::FieldList<Dimension, Tensor>    mDvDx;
  FieldSpace::FieldList<Dimension, Tensor>    mInternalDvDx;
  FieldSpace::FieldList<Dimension, Vector>    mDmassDensityDx;

  FieldSpace::FieldList<Dimension, Scalar>    mVolume;

  FieldSpace::FieldList<Dimension, std::vector<Vector> > mPairAccelerations;

  FieldSpace::FieldList<Dimension, Scalar>    mM0;
  FieldSpace::FieldList<Dimension, Vector>    mM1;
  FieldSpace::FieldList<Dimension, SymTensor> mM2;

  FieldSpace::FieldList<Dimension, Scalar>    mA0;
  FieldSpace::FieldList<Dimension, Scalar>    mA;
  FieldSpace::FieldList<Dimension, Vector>    mB;
  FieldSpace::FieldList<Dimension, Vector>    mC;
  FieldSpace::FieldList<Dimension, Tensor>    mD;
  FieldSpace::FieldList<Dimension, Vector>    mGradA0;
  FieldSpace::FieldList<Dimension, Vector>    mGradA;
  FieldSpace::FieldList<Dimension, Tensor>    mGradB;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  CSPHHydroBase();
  CSPHHydroBase(const CSPHHydroBase&);
  CSPHHydroBase& operator=(const CSPHHydroBase&);
};

}
}

#ifndef __GCCXML__
#include "CSPHHydroBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace CSPHSpace {
    template<typename Dimension> class CSPHHydroBase;
  }
}

#endif
