//---------------------------------Spheral++----------------------------------//
// SVPHFacetedHydroBase -- The SVPHFaceted hydrodynamic package for Spheral++.
//
// Created by JMO, Sun Jul 28 20:57:01 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHFacetedHydroBase_hh__
#define __Spheral_SVPHFacetedHydroBase_hh__

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
namespace SVPHSpace {

template<typename Dimension>
class SVPHFacetedHydroBase: public PhysicsSpace::GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SVPHFacetedHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                       const KernelSpace::TableKernel<Dimension>& W,
                       ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool XSVPH,
                       const PhysicsSpace::MassDensityType densityUpdate,
                       const PhysicsSpace::HEvolutionType HUpdate,
                       const Vector& xmin,
                       const Vector& xmax);

  // Destructor.
  virtual ~SVPHFacetedHydroBase();

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

  // Flag to choose whether we want to sum for density, or integrate
  // the continuity equation.
  PhysicsSpace::MassDensityType densityUpdate() const;
  void densityUpdate(const PhysicsSpace::MassDensityType type);

  // Flag to select how we want to evolve the H tensor.
  PhysicsSpace::HEvolutionType HEvolution() const;
  void HEvolution(const PhysicsSpace::HEvolutionType type);

  // Flag to determine if we're using the total energy conserving compatible energy
  // evolution scheme.
  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(const bool val);

  // Flag to determine if we're using the XSVPH algorithm.
  bool XSVPH() const;
  void XSVPH(const bool val);

  // Optionally we can provide a bounding box for use generating the mesh.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // The object defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The tessellation.
  const MeshSpace::Mesh<Dimension>& mesh() const;

  // The state field lists we're maintaining.
  // const FieldSpace::FieldList<Dimension, std::vector<Scalar> >&    A() const;
  // const FieldSpace::FieldList<Dimension, std::vector<Vector> >&    B() const;
  // const FieldSpace::FieldList<Dimension, std::vector<Tensor> >&    gradB() const;
  const FieldSpace::FieldList<Dimension, int>&       timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    volume() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    massDensitySum() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldSpace::FieldList<Dimension, Vector>&    XSVPHDeltaV() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DxDt() const;
  const FieldSpace::FieldList<Dimension, Vector>&    DvDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldSpace::FieldList<Dimension, Tensor>&    internalDvDx() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SVPHFacetedHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  // The method defining how we evolve smoothing scales.
  const NodeSpace::SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  PhysicsSpace::MassDensityType mDensityUpdate;
  PhysicsSpace::HEvolutionType mHEvolution;
  bool mCompatibleEnergyEvolution, mXSVPH;

  // Optional bounding box for generating the mesh.
  Vector mXmin, mXmax;

  // The mesh.
  typedef boost::shared_ptr<MeshSpace::Mesh<Dimension> > MeshPtr;
  MeshPtr mMeshPtr;

  // Some internal scratch fields.
  // FieldSpace::FieldList<Dimension, std::vector<Scalar> >&    mA;
  // FieldSpace::FieldList<Dimension, std::vector<Vector> >&    mB;
  // FieldSpace::FieldList<Dimension, std::vector<Tensor> >&    mGradB;
  FieldSpace::FieldList<Dimension, int>       mTimeStepMask;
  FieldSpace::FieldList<Dimension, Scalar>    mPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldSpace::FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;

  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  FieldSpace::FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldSpace::FieldList<Dimension, Scalar>    mMassDensitySum;

  FieldSpace::FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldSpace::FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldSpace::FieldList<Dimension, Vector>    mXSVPHDeltaV;

  FieldSpace::FieldList<Dimension, Vector>    mDxDt;
  FieldSpace::FieldList<Dimension, Vector>    mDvDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldSpace::FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldSpace::FieldList<Dimension, SymTensor> mDHDt;
  FieldSpace::FieldList<Dimension, Tensor>    mDvDx;
  FieldSpace::FieldList<Dimension, Tensor>    mInternalDvDx;

  FieldSpace::FieldList<Dimension, Scalar>    mVolume;

private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  SVPHFacetedHydroBase();
  SVPHFacetedHydroBase(const SVPHFacetedHydroBase&);
  SVPHFacetedHydroBase& operator=(const SVPHFacetedHydroBase&);
};

}
}

#ifndef __GCCXML__
#include "SVPHFacetedHydroBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SVPHFacetedSpace {
    template<typename Dimension> class SVPHFacetedHydroBase;
  }
}

#endif
