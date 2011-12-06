//---------------------------------Spheral++----------------------------------//
// ConstantAcceleration -- Impose a constant acceleration on a given set of
// nodes.
//
// Created by JMO, Mon Sep 27 23:01:15 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsSpace_ConstantAcceleration__
#define __Spheral_PhysicsSpace_ConstantAcceleration__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "Physics/GenericBodyForce.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace PhysicsSpace {

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
                       const NodeSpace::NodeList<Dimension>& nodeList,
                       const std::vector<int>& indicies);

  // Destructor.
  virtual ~ConstantAcceleration();

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Access the constant acceleration.
  Vector a0() const;

  // Access the NodeList.
  const NodeSpace::NodeList<Dimension>& nodeList() const;

  // Access the set of node indicies.
  const std::vector<int>& indicies() const;

private:
  //--------------------------- Public Interface ---------------------------//
  Vector ma0;
  const NodeSpace::NodeList<Dimension>* mNodeListPtr;
  const std::vector<int> mIndicies;

  // No default constructor, copying, or assignment.
  ConstantAcceleration();
  ConstantAcceleration(const ConstantAcceleration& rhs);
  ConstantAcceleration& operator=(const ConstantAcceleration& rhs);
};

}
}

#ifndef __GCCXML__
#include "ConstantAccelerationInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class ConstantAcceleration;
  }
}

#endif
