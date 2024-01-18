//---------------------------------Spheral++----------------------------------//
// LinearAcceleration -- Impose a linear acceleration as a function of position
// on a given set of nodes.
//
// Created by JMO, Wed Sep 29 15:45:23 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsSpace_LinearAcceleration__
#define __Spheral_PhysicsSpace_LinearAcceleration__

#include "Physics/GenericBodyForce.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class NodeList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class LinearAcceleration: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  LinearAcceleration(const Scalar a0,
                     const Scalar aslope);

  // Destructor.
  virtual ~LinearAcceleration();

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
  virtual std::string label() const { return "LinearAcceleration"; }

  // Access the acceleration parameters.
  Scalar a0() const;
  Scalar aslope() const;

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar ma0;
  Scalar maslope;

  // No default constructor, copying, or assignment.
  LinearAcceleration();
  LinearAcceleration(const LinearAcceleration& rhs);
  LinearAcceleration& operator=(const LinearAcceleration& rhs);
};

}

#include "LinearAccelerationInline.hh"

#endif
