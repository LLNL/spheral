//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The CRKSPH/ACRKSPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidCRKSPHHydroBase_hh__
#define __Spheral_SolidCRKSPHHydroBase_hh__

#include "CRKSPH/CRKSPHHydroBase.hh"

#include <string>

namespace Spheral {
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;
}

namespace Spheral {

template<typename Dimension>
class SolidCRKSPHHydroBase: public CRKSPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBase(DataBase<Dimension>& dataBase,
                       ArtificialViscosity<Dimension>& Q,
                       const RKOrder order,
                       const double filter,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool evolveTotalEnergy,
                       const bool XSPH,
                       const MassDensityType densityUpdate,
                       const double epsTensile,
                       const double nTensile,
                       const bool damageRelieveRubble);

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBase() = delete;
  SolidCRKSPHHydroBase(const SolidCRKSPHHydroBase&) = delete;
  SolidCRKSPHHydroBase& operator=(const SolidCRKSPHHydroBase&) = delete;

  // Destructor.
  virtual ~SolidCRKSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                            State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivatives) override;

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

  // The state field lists we're maintaining.
  const FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldList<Dimension, Scalar>&    bulkModulus() const;
  const FieldList<Dimension, Scalar>&    shearModulus() const;
  const FieldList<Dimension, Scalar>&    yieldStrength() const;
  const FieldList<Dimension, Scalar>&    plasticStrain0() const;

  // Control whether allow damaged material to have stress relieved.
  bool damageRelieveRubble() const;
  void damageRelieveRubble(bool x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidCRKSPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mDamageRelieveRubble;

private:
  //--------------------------- Private Interface ---------------------------//
  // Some internal scratch fields.
  FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldList<Dimension, Scalar> mBulkModulus;
  FieldList<Dimension, Scalar> mShearModulus;
  FieldList<Dimension, Scalar> mYieldStrength;
  FieldList<Dimension, Scalar> mPlasticStrain0;
};

}

#include "SolidCRKSPHHydroBaseInline.hh"

#endif
