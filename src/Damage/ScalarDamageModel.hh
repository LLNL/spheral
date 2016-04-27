//---------------------------------Spheral++----------------------------------//
// ScalarDamageModel -- Base class for the scalar damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendant classes.
//
// Created by JMO, Sun Oct 10 17:22:05 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ScalarDamageModel_hh__
#define __Spheral_ScalarDamageModel_hh__

#ifndef __GCCXML__
#include <vector>
#include "Geometry/GeomPlane.hh"
#include "NodeList/FluidNodeList.hh"
#else
#include "fakestl.hh"
#endif

#include "DamageModel.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class GeomPlane;
  namespace NodeSpace {
    template<typename Dimension> class FluidNodeList;
    template<typename Dimension> class SolidNodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
}

namespace Spheral {
namespace PhysicsSpace {

template<typename Dimension>
class ScalarDamageModel: 
    public DamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef GeomPlane<Dimension> Plane;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef FieldSpace::Field<Dimension, std::vector<double> > FlawStorageType;

  // Constructors, destructor.
  ScalarDamageModel(SolidMaterial::SolidNodeList<Dimension>& nodeList,
                    NodeSpace::FluidNodeList<Dimension>& damagedNodeList,
                    const double kernelExtent,
                    const double crackGrowthMultiplier,
                    const FlawStorageType& flaws);
  virtual ~ScalarDamageModel();

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
  //...........................................................................

  // Finalize method, called at the end of a time step.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBaseSpace::DataBase<Dimension>& db, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Access the damaged NodeList.
  const NodeSpace::FluidNodeList<Dimension>& damagedNodeList() const;

  // Specify whether to split or refine fully failed nodes.
  bool splitFailedNodes() const;
  void splitFailedNodes(const bool x);

  // Provide access to the state fields we maintain.
  const FieldSpace::Field<Dimension, Scalar>& strain() const;
  const FieldSpace::Field<Dimension, Scalar>& damage() const;
  const FieldSpace::Field<Dimension, Scalar>& DdamageDt() const;

  // Return a vector of undamged -> damaged node indicies.
  std::vector<int> undamagedToDamagedNodeIndicies() const;

  // The set of boundary planes to respect if creating new nodes.
  std::vector<Plane>& boundPlanes();

  // Determine if the given position violates any of the specified boundary
  // planes.
  bool positionOutOfBounds(const Vector& r) const;

  //**************************************************************************
  // Restart methods.
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //**************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  NodeSpace::FluidNodeList<Dimension>* mDamagedNodeListPtr;
  bool mSplitFailedNodes;

  FieldSpace::Field<Dimension, Scalar> mStrain;
  FieldSpace::Field<Dimension, Scalar> mDamage;
  FieldSpace::Field<Dimension, Scalar> mDdamageDt;

  FieldSpace::Field<Dimension, Scalar> mMass0;
  FieldSpace::Field<Dimension, int> mUndamagedToDamagedIndex;

  std::vector<Plane> mBoundPlanes;
#endif

  // No default constructor, copying or assignment.
  ScalarDamageModel();
  ScalarDamageModel(const ScalarDamageModel&);
  ScalarDamageModel& operator=(const ScalarDamageModel&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class ScalarDamageModel;
  }
}

#endif

