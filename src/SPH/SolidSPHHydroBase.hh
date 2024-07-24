//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPHHydroBase_hh__
#define __Spheral_SolidSPHHydroBase_hh__

#include <float.h>
#include <string>

#include "SPH/SPHHydroBase.hh"

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
class SolidSPHHydroBase: public SPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                    DataBase<Dimension>& dataBase,
                    ArtificialViscosity<Dimension>& Q,
                    const TableKernel<Dimension>& W,
                    const TableKernel<Dimension>& WPi,
                    const TableKernel<Dimension>& WGrad,
                    const double filter,
                    const double cfl,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool gradhCorrection,
                    const bool XSPH,
                    const bool correctVelocityGradient,
                    const bool sumMassDensityOverAllNodeLists,
                    const MassDensityType densityUpdate,
                    const HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const bool damageRelieveRubble,
                    const bool strengthInDamage,
                    const Vector& xmin,
                    const Vector& xmax);

  // Destructor.
  virtual ~SolidSPHHydroBase();

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

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Gradient kernel
  const TableKernel<Dimension>& GradKernel() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldList<Dimension, Scalar>& bulkModulus() const;
  const FieldList<Dimension, Scalar>& shearModulus() const;
  const FieldList<Dimension, Scalar>& yieldStrength() const;
  const FieldList<Dimension, Scalar>& plasticStrain0() const;
  const FieldList<Dimension, SymTensor>& Hfield0() const;

  // Control whether allow damaged material to have stress relieved.
  bool damageRelieveRubble() const;
  void damageRelieveRubble(bool x);

  // Do we allow damaged material to have strength?
  bool strengthInDamage() const;
  void strengthInDamage(bool x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidSPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mDamageRelieveRubble, mStrengthInDamage;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mGradKernel;   // Gradient kernel

  // Some internal scratch fields.
  FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldList<Dimension, Scalar> mBulkModulus;
  FieldList<Dimension, Scalar> mShearModulus;
  FieldList<Dimension, Scalar> mYieldStrength;
  FieldList<Dimension, Scalar> mPlasticStrain0;
  FieldList<Dimension, SymTensor> mHfield0;

  // No default constructor, copying, or assignment.
  SolidSPHHydroBase();
  SolidSPHHydroBase(const SolidSPHHydroBase&);
  SolidSPHHydroBase& operator=(const SolidSPHHydroBase&);
};

}

#include "SolidSPHHydroBaseInline.hh"

#endif
