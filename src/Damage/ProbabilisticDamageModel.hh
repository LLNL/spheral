//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamageModel
//
// A damage model based on Weibull statistics that uses volume based
// probabilities per node to decide when damage starts to accrue.  Should
// generate similar results to the classic Benz-Asphaug (Grady-Kipp) model
// without generating explicit flaws.  Also appropriate for use with varying
// resolution materials.
//
// Created by JMO, Tue Apr 13 15:58:08 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_ProbabilisticDamageModel_hh__
#define __Spheral_ProbabilisticDamageModel_hh__

#include "TensorDamageModel.hh"  // For now, so we pick up the enums

template<typename Dimension>
class ProbabilisticDamageModel: public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using FlawStorageType = Field<Dimension, std::vector<double> >;

  // Constructors, destructor.
  ProbabilisticDamageModel(SolidNodeList<Dimension>& nodeList,
                           const TableKernel<Dimension>& W,
                           const double crackGrowthMultiplier,
                           const DamageCouplingAlgorithm damageCouplingAlgorithm,
                           const TensorStrainAlgorithm strainAlgorithm,
                           const bool damageInCompression);
  virtual ~ProbabilisticDamageModel();

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

  //............................................................................
  // Accessors for state
  TensorStrainAlgorithm strainAlgorithm() const;
  bool damageInCompression() const;
  const Field<Dimension, unsigned>& numFlaws() const;
  const Field<Dimension, unsigned>& numFlawsActivated() const;
  const Field<Dimension, Scalar>& currentFlaw() const;
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;
  const Field<Dimension, SymTensor>& strain() const;
  const Field<Dimension, Scalar>& DdamageDt() const;

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainAlgorithm;
  bool mDamageInCompression;
  Field<Dimension, unsigned> mNumFlaws, mNumFlawsActivated;
  Field<Dimension, Scalar> mCurrentFlaw, mYoungsModulus, mLongitudinalSoundSpeed, mDdamageDt;
  Field<Dimension, SymTensor> mStrain, mEffectiveStrain;

  // No default constructor, copying or assignment.
  ProbabilisticDamageModel();
  ProbabilisticDamageModel(const ProbabilisticDamageModel&);
  ProbabilisticDamageModel& operator=(const ProbabilisticDamageModel&);
};

}

#include "ProbabilisticDamageModelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ProbabilisticDamageModel;
}

#endif

