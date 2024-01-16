//---------------------------------Spheral++----------------------------------//
// ConstantAcceleration -- Impose a constant acceleration on a given set of
// nodes.
//
// Created by JMO, Mon Sep 27 23:01:15 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsSpace_ConstantAcceleration__
#define __Spheral_PhysicsSpace_ConstantAcceleration__

#include "Physics/GenericBodyForce.hh"
#include "Field/Field.hh"

#include <vector>
#include <memory>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class NodeList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class ConstantAcceleration: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  ConstantAcceleration(const Vector a0,
                       const NodeList<Dimension>& nodeList,
                       const std::vector<int>& indices);
  ConstantAcceleration(const Vector a0,
                       const NodeList<Dimension>& nodeList);

  // Destructor.
  virtual ~ConstantAcceleration();

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Required label for Physics interface.
  virtual std::string label() const { return "ConstantAcceleration"; }

  // Access the constant acceleration.
  Vector a0() const;

  // Access the NodeList.
  const NodeList<Dimension>& nodeList() const;

  // Access the set of node flags.
  const Field<Dimension, int>& flags() const;

private:
  //--------------------------- Public Interface ---------------------------//
  Vector ma0;
  const NodeList<Dimension>* mNodeListPtr;
  std::shared_ptr<Field<Dimension, int>> mFlagsPtr;

  // No default constructor, copying, or assignment.
  ConstantAcceleration();
  ConstantAcceleration(const ConstantAcceleration& rhs);
  ConstantAcceleration& operator=(const ConstantAcceleration& rhs);
};

}

#include "ConstantAccelerationInline.hh"

#endif
