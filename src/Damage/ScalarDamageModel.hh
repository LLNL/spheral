//---------------------------------Spheral++----------------------------------//
// ScalarDamageModel -- Base class for the scalar damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendant classes.
//
// Created by JMO, Sun Oct 10 17:22:05 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ScalarDamageModel_hh__
#define __Spheral_ScalarDamageModel_hh__

#include "Geometry/GeomPlane.hh"
#include "NodeList/FluidNodeList.hh"
#include "DamageModel.hh"

#include <vector>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class GeomPlane;
  template<typename Dimension> class FluidNodeList;
  template<typename Dimension> class SolidNodeList;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
}

namespace Spheral {

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
  typedef Field<Dimension, std::vector<double> > FlawStorageType;

  // Constructors, destructor.
  ScalarDamageModel(SolidNodeList<Dimension>& nodeList,
                    FluidNodeList<Dimension>& damagedNodeList,
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
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);
  //...........................................................................

  // Finalize method, called at the end of a time step.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& db, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Access the damaged NodeList.
  const FluidNodeList<Dimension>& damagedNodeList() const;

  // Specify whether to split or refine fully failed nodes.
  bool splitFailedNodes() const;
  void splitFailedNodes(const bool x);

  // Provide access to the state fields we maintain.
  const Field<Dimension, Scalar>& strain() const;
  const Field<Dimension, Scalar>& damage() const;
  const Field<Dimension, Scalar>& DdamageDt() const;

  // Return a vector of undamged -> damaged node indicies.
  std::vector<int> undamagedToDamagedNodeIndicies() const;

  // The set of boundary planes to respect if creating new nodes.
  std::vector<Plane>& boundPlanes();

  // Determine if the given position violates any of the specified boundary
  // planes.
  bool positionOutOfBounds(const Vector& r) const;

  //**************************************************************************
  // Restart methods.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //**************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  FluidNodeList<Dimension>* mDamagedNodeListPtr;
  bool mSplitFailedNodes;

  Field<Dimension, Scalar> mStrain;
  Field<Dimension, Scalar> mDamage;
  Field<Dimension, Scalar> mDdamageDt;

  Field<Dimension, Scalar> mMass0;
  Field<Dimension, int> mUndamagedToDamagedIndex;

  std::vector<Plane> mBoundPlanes;

  // No default constructor, copying or assignment.
  ScalarDamageModel();
  ScalarDamageModel(const ScalarDamageModel&);
  ScalarDamageModel& operator=(const ScalarDamageModel&);
};

}

#endif

