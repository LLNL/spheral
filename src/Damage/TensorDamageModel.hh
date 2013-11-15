//---------------------------------Spheral++----------------------------------//
// TensorDamageModel -- Base class for the tensor damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendant classes.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Thu Sep 29 15:42:05 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorDamageModel_hh__
#define __Spheral_TensorDamageModel_hh__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "DamageModel.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class SolidNodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace PhysicsSpace {

// Enum for selecting the method of defining the tensor strain.
enum TensorStrainAlgorithm {
  BenzAsphaug = 0,
  StrainHistory = 1,
  MeloshRyanAsphaug = 2,
  PlasticStrain = 3,
  PseudoPlasticStrain = 4,
};

// Enum for selecting the method of defining the effective tensor damage.
enum EffectiveDamageAlgorithm {
  Copy = 0,
  Max = 1,
  Sampled = 2,
};

template<typename Dimension>
class TensorDamageModel: 
    public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef FieldSpace::Field<Dimension, std::vector<double> > FlawStorageType;

  // Constructors, destructor.
  TensorDamageModel(SolidMaterial::SolidNodeList<Dimension>& nodeList,
                    const TensorStrainAlgorithm strainAlgorithm,
                    const EffectiveDamageAlgorithm effDamageAlgorithm,
                    const bool useDamageGradient,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const double crackGrowthMultiplier,
                    const EffectiveFlawAlgorithm flawAlgorithm,
                    const double criticalDamageThreshold,
                    const FlawStorageType& flaws);
  virtual ~TensorDamageModel();

  //...........................................................................
  // Provide the required physics package interface.
  // Compute the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register our state.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);
  //...........................................................................

  // Provide access to the state fields we maintain.
  const FieldSpace::Field<Dimension, SymTensor>& strain() const;
  const FieldSpace::Field<Dimension, SymTensor>& effectiveStrain() const;
  const FieldSpace::Field<Dimension, Scalar>& DdamageDt() const;
  const FieldSpace::Field<Dimension, SymTensor>& newEffectiveDamage() const;
  const FieldSpace::Field<Dimension, Vector>& newDamageGradient() const;

  // The algorihms being used to update the strain and effective damage.
  TensorStrainAlgorithm strainAlgorithm() const;
  EffectiveDamageAlgorithm effectiveDamageAlgorithm() const;

  // Flag to determine if we compute the gradient of the damage at the start 
  // of a timestep.
  bool useDamageGradient() const;
  void useDamageGradient(const bool x);

  // The critical damage threshold for not setting the time step.
  double criticalDamageThreshold() const; 
  void criticalDamageThreshold(const double x);

  //**************************************************************************
  // Restart methods.
  virtual std::string label() const { return "TensorDamageModel"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //**************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
#ifndef __GCCXML__
  FieldSpace::Field<Dimension, SymTensor> mStrain;
  FieldSpace::Field<Dimension, SymTensor> mEffectiveStrain;
  FieldSpace::Field<Dimension, Scalar> mDdamageDt;
  FieldSpace::Field<Dimension, SymTensor> mNewEffectiveDamage;
  FieldSpace::Field<Dimension, Vector> mNewDamageGradient;
#endif

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainAlgorithm;
  EffectiveDamageAlgorithm mEffDamageAlgorithm;
  double mCriticalDamageThreshold;
  bool mUseDamageGradient;

  // No default constructor, copying or assignment.
  TensorDamageModel();
  TensorDamageModel(const TensorDamageModel&);
  TensorDamageModel& operator=(const TensorDamageModel&);
};

}
}

#ifndef __GCCXML__
#include "TensorDamageModelInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class TensorDamageModel;
  }
}

#endif

