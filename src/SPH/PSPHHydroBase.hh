//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_PSPHHydroBase_hh__
#define __Spheral_PSPHHydroBase_hh__

#include <string>

#include "SPHHydroBase.hh"

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
class PSPHHydroBase: public SPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  PSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                DataBase<Dimension>& dataBase,
                ArtificialViscosity<Dimension>& Q,
                const TableKernel<Dimension>& W,
                const TableKernel<Dimension>& WPi,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const bool correctVelocityGradient,
                const bool HopkinsConductivity,
                const bool sumMassDensityOverAllNodeLists,
                const MassDensityType densityUpdate,
                const HEvolutionType HUpdate,
                const Vector& xmin,
                const Vector& xmax);

  // Destructor.
  virtual ~PSPHHydroBase();

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Pre-step initializations.
  virtual 
  void preStepInitialize(const DataBase<Dimension>& dataBase, 
                         State<Dimension>& state,
                         StateDerivatives<Dimension>& derivatives) override;

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

  // Post-state update. For PSPH this is where we recompute the PSPH pressure and corrections.
  virtual 
  void postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Flag determining if we're applying Hopkins 2014 conductivity.
  bool HopkinsConductivity() const;
  void HopkinsConductivity(bool val);

  // The state field lists we're maintaining.
  const FieldList<Dimension, Scalar>&    gamma() const;
  const FieldList<Dimension, Scalar>&    PSPHcorrection() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "PSPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  bool mHopkinsConductivity;

  //PSPH Fields
  FieldList<Dimension, Scalar>    mGamma;
  FieldList<Dimension, Scalar>    mPSPHcorrection;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  PSPHHydroBase();
  PSPHHydroBase(const PSPHHydroBase&);
  PSPHHydroBase& operator=(const PSPHHydroBase&);
};

}

#include "PSPHHydroBaseInline.hh"

#endif
