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

#include "DamageModel.hh"
#include "TensorDamageModel.hh"             // For now, so we pick up the enums
#include "Utilities/uniform_random.hh"

namespace Spheral {

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

  // Constructors, destructor.
  ProbabilisticDamageModel(SolidNodeList<Dimension>& nodeList,
                           const TableKernel<Dimension>& W,
                           const double kWeibull,
                           const double mWeibull,
                           const size_t seed,
                           const size_t minFlawsPerNode,
                           const double crackGrowthMultiplier,
                           const double volumeMultiplier,
                           const DamageCouplingAlgorithm damageCouplingAlgorithm,
                           const TensorStrainAlgorithm strainAlgorithm,
                           const bool damageInCompression);
  virtual ~ProbabilisticDamageModel();

  //............................................................................
  // Override the Physics package interface.

  // Initialize once when the problem is starting up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

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

  //............................................................................
  // Accessors for state
  TensorStrainAlgorithm strainAlgorithm() const;
  bool damageInCompression() const;
  double kWeibull() const;
  double mWeibull() const;
  double volumeMultiplier() const;
  double Vmin() const;
  double Vmax() const;
  size_t seed() const;
  size_t minFlawsPerNode() const;
  const Field<Dimension, unsigned>& numFlaws() const;
  const Field<Dimension, Scalar>& minFlaw() const;
  const Field<Dimension, Scalar>& maxFlaw() const;
  const Field<Dimension, Scalar>& initialVolume() const;
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;
  const Field<Dimension, SymTensor>& strain() const;
  const Field<Dimension, SymTensor>& effectiveStrain() const;
  const Field<Dimension, Scalar>& DdamageDt() const;
  const Field<Dimension, uniform_random>& randomGenerator() const;

  // Optionally the user can provide a mask to prevent damage modeling on some points.
  const Field<Dimension, int>& mask() const;
  void mask(const Field<Dimension, int>& val);

  //............................................................................
  // Restart methods.
  virtual std::string label() const { return "ProbabilisticDamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainAlgorithm;
  bool mDamageInCompression;
  double mkWeibull, mmWeibull, mVolumeMultiplier, mVmin, mVmax;
  size_t mSeed, mMinFlawsPerNode;
  Field<Dimension, unsigned> mNumFlaws;
  Field<Dimension, Scalar> mMinFlaw, mMaxFlaw, mInitialVolume, mYoungsModulus, mLongitudinalSoundSpeed, mDdamageDt;
  Field<Dimension, SymTensor> mStrain, mEffectiveStrain;
  Field<Dimension, uniform_random> mRandomGenerator;
  Field<Dimension, int> mMask;

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

