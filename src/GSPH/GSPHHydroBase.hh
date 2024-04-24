//---------------------------------Spheral++----------------------------------//
// GSPHHydroBase -- The Godunov SPH hydrodynamic package for Spheral++.
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_GSPHHydroBase_hh__
#define __Spheral_GSPHHydroBase_hh__

#include <string>

#include "GSPH/GenericRiemannHydro.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class RiemannSolverBase;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class GSPHHydroBase: public GenericRiemannHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;

  using TimeStepType = typename GenericRiemannHydro<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename GenericRiemannHydro<Dimension>::ConstBoundaryIterator;

  // Constructors.
  GSPHHydroBase(DataBase<Dimension>& dataBase,
                RiemannSolverBase<Dimension>& riemannSolver,
                const TableKernel<Dimension>& W,
                const Scalar epsDiffusionCoeff,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const bool correctVelocityGradient,
                const GradientType gradType,
                const MassDensityType densityUpdate,
                const double epsTensile,
                const double nTensile,
                const Vector& xmin,
                const Vector& xmax);

  
  // No default constructor, copying, or assignment.
  GSPHHydroBase() = delete;
  GSPHHydroBase(const GSPHHydroBase&) = delete;
  GSPHHydroBase& operator=(const GSPHHydroBase&) = delete;

  // Destructor.
  virtual ~GSPHHydroBase();

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
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

  void
  computeMCorrection(const typename Dimension::Scalar time,
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


  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "GSPHHydroBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************
  
private:
  FieldList<Dimension, Scalar>    mDmassDensityDt;
};

}

#include "GSPHHydroBaseInline.hh"

#endif
