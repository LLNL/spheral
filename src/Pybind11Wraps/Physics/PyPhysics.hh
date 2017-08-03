//------------------------------------------------------------------------------
// Trampoline classes for abstract Physics interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyPhysics__
#define __Spheral_PyPhysics__

#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// PyPhysics
//------------------------------------------------------------------------------
template<typename Dimension, class PhysicsBase>
class PyPhysics: public PhysicsBase {
public:
  using PhysicsBase::PhysicsBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivs) const override {
    PYBIND11_OVERLOAD(void,                // Return type
                           PhysicsBase,         // Parent class
                           evaluateDerivatives, // name of method
                           time, dt, dataBase, state, derivs   // arguments
                           );
  }

  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override {
    PYBIND11_OVERLOAD(TimeStepType,        // Return type
                           PhysicsBase,         // Parent class
                           dt,                  // name of method
                           dataBase, state, derivs, currentTime   // arguments
                           );
  }

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override {
    PYBIND11_OVERLOAD(void,                // Return type
                           PhysicsBase,         // Parent class
                           registerState,       // name of method
                           dataBase, state   // arguments
                           );
  }

  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,                // Return type
                           PhysicsBase,         // Parent class
                           registerDerivatives, // name of method
                           dataBase, derivs   // arguments
                           );
  }

  virtual std::string label() const override {
    PYBIND11_OVERLOAD(std::string,         // Return type
                           PhysicsBase,         // Parent class
                           label // name of method
                           );
  }

  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      applyGhostBoundaries, // name of method
                      state, derivs        // arguments
                      );
  }

  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      enforceBoundaries, // name of method
                      state, derivs        // arguments
                      );
  }

  virtual void initialize(const Scalar time,
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase,
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      initialize, // name of method
                      time, dt, dataBase, state, derivs        // arguments
                      );
  }

  virtual void finalize(const Scalar time,
                        const Scalar dt,
                        DataBase<Dimension>& dataBase,
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      finalize, // name of method
                      time, dt, dataBase, state, derivs        // arguments
                      );
  }

  virtual void finalizeDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivs) const override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      finalizeDerivatives, // name of method
                      time, dt, dataBase, state, derivs        // arguments
                      );
  }

  virtual void postStateUpdate(const DataBase<Dimension>& dataBase,
                               State<Dimension>& state,
                               const StateDerivatives<Dimension>& derivs) const override {
    PYBIND11_OVERLOAD(void,                // Return type
                      PhysicsBase,         // Parent class
                      postStateUpdate, // name of method
                      dataBase, state, derivs        // arguments
                      );
  }

  virtual bool requireConnectivity() const override {
    PYBIND11_OVERLOAD(bool,                // Return type
                      PhysicsBase,         // Parent class
                      requireConnectivity, // name of method
                      );
  }

  virtual bool requireGhostConnectivity() const override {
    PYBIND11_OVERLOAD(bool,                // Return type
                      PhysicsBase,         // Parent class
                      requireGhostConnectivity, // name of method
                      );
  }

  virtual Scalar extraEnergy() const override {
    PYBIND11_OVERLOAD(Scalar,                // Return type
                      PhysicsBase,         // Parent class
                      extraEnergy, // name of method
                      );
  }

  virtual Vector extraMomentum() const override {
    PYBIND11_OVERLOAD(Vector,                // Return type
                      PhysicsBase,         // Parent class
                      extraMomentum, // name of method
                      );
  }

};

}
}

#endif
