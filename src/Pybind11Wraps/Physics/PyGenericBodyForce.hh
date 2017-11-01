//------------------------------------------------------------------------------
// Trampoline classes for GenericBodyForce interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyGenericBodyForce__
#define __Spheral_PyGenericBodyForce__

#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Physics/GenericBodyForce.hh"
#include "PyAbstractPhysics.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// PyGenericBodyForce
//------------------------------------------------------------------------------
template<typename Dimension>
class PyGenericBodyForce: public PyAbstractPhysics<Dimension, GenericBodyForce<Dimension>> {
public:
  typedef GenericBodyForce<Dimension> PhysicsBase;
  typedef PyAbstractPhysics<Dimension, PhysicsBase> PyAP;
  using PyAP::PyAP;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

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

};

}
}

#endif
