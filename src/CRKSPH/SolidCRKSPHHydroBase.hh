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
template<typename Dimension> class SmoothingScaleBase;
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
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                       DataBase<Dimension>& dataBase,
                       ArtificialViscosity<Dimension>& Q,
                       const RKOrder order,
                       const double filter,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool evolveTotalEnergy,
                       const bool XSPH,
                       const MassDensityType densityUpdate,
                       const HEvolutionType HUpdate,
                       const double epsTensile,
                       const double nTensile,
                       const bool damageRelieveRubble);

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
  const FieldList<Dimension, SymTensor>& Hfield0() const;

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
  FieldList<Dimension, SymTensor> mHfield0;

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBase();
  SolidCRKSPHHydroBase(const SolidCRKSPHHydroBase&);
  SolidCRKSPHHydroBase& operator=(const SolidCRKSPHHydroBase&);
};

}

#include "SolidCRKSPHHydroBaseInline.hh"

#endif
