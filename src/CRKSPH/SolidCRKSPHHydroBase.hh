//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
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

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                       ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                       const KernelSpace::TableKernel<Dimension>& W,
                       const KernelSpace::TableKernel<Dimension>& WPi,
                       const double filter,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool evolveTotalEnergy,
                       const bool XSPH,
                       const PhysicsSpace::MassDensityType densityUpdate,
                       const PhysicsSpace::HEvolutionType HUpdate,
                       const CRKSPHSpace::CRKOrder correctionOrder,
                       const CRKSPHSpace::CRKVolumeType volumeType,
                       const double epsTensile,
                       const double nTensile,
                       const bool damageRelieveRubble);

  // Destructor.
  virtual ~SolidCRKSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
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
  const FieldSpace::FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>& bulkModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& shearModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& yieldStrength() const;
  const FieldSpace::FieldList<Dimension, Scalar>& plasticStrain0() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& Hfield0() const;
  const FieldSpace::FieldList<Dimension, int>& fragIDs() const;

  // Control whether allow damaged material to have stress relieved.
  bool damageRelieveRubble() const;
  void damageRelieveRubble(const bool x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SolidCRKSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  bool mDamageRelieveRubble;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldSpace::FieldList<Dimension, Scalar> mBulkModulus;
  FieldSpace::FieldList<Dimension, Scalar> mShearModulus;
  FieldSpace::FieldList<Dimension, Scalar> mYieldStrength;
  FieldSpace::FieldList<Dimension, Scalar> mPlasticStrain0;
  FieldSpace::FieldList<Dimension, SymTensor> mHfield0;
  FieldSpace::FieldList<Dimension, int> mFragIDs;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBase();
  SolidCRKSPHHydroBase(const SolidCRKSPHHydroBase&);
  SolidCRKSPHHydroBase& operator=(const SolidCRKSPHHydroBase&);
};

}

#include "SolidCRKSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SolidCRKSPHHydroBase;
}

#endif
