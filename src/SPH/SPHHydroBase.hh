//---------------------------------Spheral++----------------------------------//
// SPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHHydroBase_hh__
#define __Spheral_SPHHydroBase_hh__

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
namespace SPHSpace {

template<typename Dimension>
class SPHHydroBase: public PhysicsSpace::GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
               const KernelSpace::TableKernel<Dimension>& W,
               const KernelSpace::TableKernel<Dimension>& WPi,
               ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool gradhCorrection,
               const bool XSPH,
               const bool correctVelocityGradient,
               const bool sumMassDensityOverAllNodeLists,
               const PhysicsSpace::MassDensityType densityUpdate,
               const PhysicsSpace::HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax);

  // Destructor.
  virtual ~SPHHydroBase();

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

  // A method to fill in the volume in the State, optionally enforcing
  // boundary conditions.
  void updateVolume(State<Dimension>& state,
                    const bool boundaries) const;

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

  // Flag to determine if we're applying the linear correction for the velocity gradient.
  bool correctVelocityGradient() const;
  void correctVelocityGradient(const bool val);

  // Flag to determine if the sum density definition extends over neighbor NodeLists.
  bool sumMassDensityOverAllNodeLists() const;
  void sumMassDensityOverAllNodeLists(const bool val);

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile() const;
  void epsilonTensile(const Scalar val);

  Scalar nTensile() const;
  void nTensile(const Scalar val);

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // The object defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, int>&       timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    nodeScale() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    volume() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    omegaGradh() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    massDensitySum() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    normalization() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldSpace::FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    M() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    localM() const;
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
  virtual std::string label() const { return "SPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, std::string pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  PhysicsSpace::MassDensityType mDensityUpdate;
  PhysicsSpace::HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mGradhCorrection, mXSPH, mCorrectVelocityGradient, mSumMassDensityOverAllNodeLists;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, int>       mTimeStepMask;
  FieldSpace::FieldList<Dimension, Scalar>    mNodeScale;

  FieldSpace::FieldList<Dimension, Scalar>    mPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldSpace::FieldList<Dimension, Scalar>    mOmegaGradh;
  FieldSpace::FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  FieldSpace::FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mMassDensitySum;
  FieldSpace::FieldList<Dimension, Scalar>    mNormalization;

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
  FieldSpace::FieldList<Dimension, Tensor>    mM;
  FieldSpace::FieldList<Dimension, Tensor>    mLocalM;

  FieldSpace::FieldList<Dimension, Scalar>    mVolume;

  FieldSpace::FieldList<Dimension, std::vector<Vector> > mPairAccelerations;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  SPHHydroBase();
  SPHHydroBase(const SPHHydroBase&);
  SPHHydroBase& operator=(const SPHHydroBase&);
};

}
}

#ifndef __GCCXML__
#include "SPHHydroBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SPHSpace {
    template<typename Dimension> class SPHHydroBase;
  }
}

#endif
