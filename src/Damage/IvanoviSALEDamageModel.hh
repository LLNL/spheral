//---------------------------------Spheral++----------------------------------//
// IvanoviSALEDamageModel
//
// The Ivanov damage model, hopefully close to how it's implemented in iSALE.
// This damage model is most appropriate for rocky materials.
//
// Refs:
// 
// Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
//   Meteoritics & Planetary Science, 39(2), 217. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x
//
// Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
//   internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040
//
// Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
//   Mining Sciences & Geomechanics Abstracts, 4(3):269
//
// Created by JMO, Sat Jun 26 11:35:44 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_IvanoviSALEDamageModel_hh__
#define __Spheral_IvanoviSALEDamageModel_hh__

#include "DamageModel.hh"
#include "TensorDamageModel.hh"             // For now, so we pick up the enums

namespace Spheral {

template<typename Dimension>
class IvanoviSALEDamageModel: public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors, destructor.
  IvanoviSALEDamageModel(SolidNodeList<Dimension>& nodeList,
                         const TableKernel<Dimension>& W,
                         const double minPlasticFailure,
                         const double plasticFailurePressureSlope,
                         const double plasticFailurePressureOffset,
                         const double tensileFailureStress,
                         const double crackGrowthMultiplier,
                         const DamageCouplingAlgorithm damageCouplingAlgorithm,
                         const double criticalDamageThreshold,
                         const Field<Dimension, int>& mask);
  virtual ~IvanoviSALEDamageModel();

  //............................................................................
  // Override the Physics package interface.
  // Compute the derivatives.
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;
  //............................................................................
  // Accessors for state
  double minPlasticFailure() const;
  double plasticFailurePressureSlope() const;
  double plasticFailurePressureOffset() const;
  double tensileFailureStress() const;
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;
  const Field<Dimension, SymTensor>& strain() const;
  const Field<Dimension, SymTensor>& effectiveStrain() const;
  const Field<Dimension, Scalar>& DdamageDt() const;

  // Optionally the user can provide a mask to prevent damage modeling on some points.
  const Field<Dimension, int>& mask() const;
  void mask(const Field<Dimension, int>& val);

  // Optionally ignore timestep votes for material beyond a damage threshold
  double criticalDamageThreshold() const;
  void criticalDamageThreshold(const double val);

  //............................................................................
  // Restart methods.
  virtual std::string label() const override { return "IvanoviSALEDamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;

private:
  //--------------------------- Private Interface ---------------------------//
  double mEpsPfb, mB, mPc, mTensileFailureStress, mCriticalDamageThreshold;
  Field<Dimension, int> mMask;
  Field<Dimension, Scalar> mYoungsModulus, mLongitudinalSoundSpeed, mDdamageDt;
  Field<Dimension, SymTensor> mStrain, mEffectiveStrain;

  // No default constructor, copying or assignment.
  IvanoviSALEDamageModel();
  IvanoviSALEDamageModel(const IvanoviSALEDamageModel&);
  IvanoviSALEDamageModel& operator=(const IvanoviSALEDamageModel&);
};

}

#include "IvanoviSALEDamageModelInline.hh"

#endif

