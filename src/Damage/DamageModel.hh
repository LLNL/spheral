//---------------------------------Spheral++----------------------------------//
// DamageModel -- Base class for the damage physics models.
// This class just provides the basic interface for damage models, and does 
// not fill out the complete physics package interface.
//
// Created by JMO, Thu Sep 29 13:31:57 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamageModel_hh__
#define __Spheral_DamageModel_hh__

#include "Physics/Physics.hh"
#include "DataOutput/registerWithRestart.hh"

#include <vector>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class SolidNodeList;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  template<typename Dimension> class TableKernel;
  class FileIO;
}

namespace Spheral {

enum class EffectiveFlawAlgorithm {
  FullSpectrumFlaws = 0,
  MinFlaw = 1,
  MaxFlaw = 2,
  InverseSumFlaws = 3,
  SampledFlaws = 4,
};

template<typename Dimension>
class DamageModel: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef Field<Dimension, std::vector<double> > FlawStorageType;

  // Constructors, destructor.
  DamageModel(SolidNodeList<Dimension>& nodeList,
              const TableKernel<Dimension>& W,
              const double crackGrowthMultiplier,
              const EffectiveFlawAlgorithm flawAlgorithm,
              const FlawStorageType& flaws);
  virtual ~DamageModel();

  // Compute the generic Grady-Kipp (ala Benz-Asphaug) scalar damage time 
  // derivative.
  virtual 
  void computeScalarDDDt(const DataBase<Dimension>& dataBase,
                         const State<Dimension>& state,
                         const Scalar time,
                         const Scalar dt,
                         Field<Dimension, Scalar>& DDDt) const;

  //...........................................................................
  // Provide a subset of the required physics package interface.
  // Descendant classes must complete the set!
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  virtual 
  void postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  virtual
  TimeStepType dt(const DataBase<Dimension>& dataBase, 
                  const State<Dimension>& state,
                  const StateDerivatives<Dimension>& derivs,
                  const Scalar currentTime) const = 0;
  //...........................................................................

  // Optional method to cull the set of flaws to the single weakest one on
  // each point.
  void cullToWeakestFlaws();

  // Get the set of flaw activation energies for the given node index.
  const std::vector<double> flawsForNode(const size_t index) const;

  // Access the SolidNodeList we're damaging.
  SolidNodeList<Dimension>& nodeList();
  const SolidNodeList<Dimension>& nodeList() const;

  // Access the kernel.
  const TableKernel<Dimension>& kernel() const;

  // Important local parameters.
  double crackGrowthMultiplier() const;
  EffectiveFlawAlgorithm effectiveFlawAlgorithm() const;

  // Allow the user to specify a set of nodes to be excluded from damage.
  std::vector<int> excludeNodes() const;
  void excludeNodes(std::vector<int> ids);

  // Provide access to the state fields we maintain.
  const Field<Dimension, Scalar>& youngsModulus() const;
  const Field<Dimension, Scalar>& longitudinalSoundSpeed() const;

  // Access the flaw field.
  const FlawStorageType& flaws() const;
  FlawStorageType& flaws();

  // The effective flaw for each node.
  const Field<Dimension, Scalar>& effectiveFlaws() const;

  // Compute a Field with the sum of the activation energies per node.
  Field<Dimension, Scalar> sumActivationEnergiesPerNode() const;

  // Compute a Field with the number of flaws per node.
  Field<Dimension, Scalar> numFlawsPerNode() const;

  // The effective critical number of nodes per smoothing scale, below which we
  // assume all flaws are active on a node.
  double criticalNodesPerSmoothingScale() const;
  void criticalNodesPerSmoothingScale(double x);

  //**************************************************************************
  // Restart methods.
  virtual std::string label() const override { return "DamageModel"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //**************************************************************************

protected:
  //-------------------------- Protected Interface --------------------------//
  FlawStorageType mFlaws;
  Field<Dimension, Scalar> mEffectiveFlaws;

private:
  //--------------------------- Private Interface ---------------------------//
  SolidNodeList<Dimension>& mNodeList;
  const TableKernel<Dimension>& mW;
  double mCrackGrowthMultiplier;
  EffectiveFlawAlgorithm mEffectiveFlawAlgorithm;

  Field<Dimension, Scalar> mYoungsModulus;
  Field<Dimension, Scalar> mLongitudinalSoundSpeed;

  Field<Dimension, int> mExcludeNode;

  double mCriticalNodesPerSmoothingScale;

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying or assignment.
  DamageModel();
  DamageModel(const DamageModel&);
  DamageModel& operator=(const DamageModel&);
};

}

#include "DamageModelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DamageModel;
}

#endif

