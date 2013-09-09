//---------------------------------Spheral++----------------------------------//
// CompatibleFaceSpecificThermalEnergyPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CompatibleFaceSpecificThermalEnergyPolicy_hh__
#define __Spheral_CompatibleFaceSpecificThermalEnergyPolicy_hh__

#include <string>

#include "DataBase/IncrementState.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class Field;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}

template<typename Dimension>
class CompatibleFaceSpecificThermalEnergyPolicy: 
    public IncrementState<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename FieldUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;
  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors, destructor.
  CompatibleFaceSpecificThermalEnergyPolicy(const KernelSpace::TableKernel<Dimension>& W,
                                            const DataBaseSpace::DataBase<Dimension>& dataBase,
                                            ConstBoundaryIterator boundaryBegin,
                                            ConstBoundaryIterator boundaryEnd);
  virtual ~CompatibleFaceSpecificThermalEnergyPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // If the derivative stored values for the pair-accelerations has not been updated,
  // we need to just time advance normally.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    IncrementState<Dimension, Scalar>::update(key,
                                              state,
                                              derivs,
                                              multiplier,
                                              t,
                                              dt);
  }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const KernelSpace::TableKernel<Dimension>& mW;
  const DataBaseSpace::DataBase<Dimension>& mDataBase;
  ConstBoundaryIterator mBoundaryBegin, mBoundaryEnd;
  static bool mFired;

  CompatibleFaceSpecificThermalEnergyPolicy(const CompatibleFaceSpecificThermalEnergyPolicy& rhs);
  CompatibleFaceSpecificThermalEnergyPolicy& operator=(const CompatibleFaceSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CompatibleFaceSpecificThermalEnergyPolicy;
}

#endif
