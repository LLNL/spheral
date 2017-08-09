//------------------------------------------------------------------------------
// Trampoline classes for abstract Physics interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyGenericHydro__
#define __Spheral_PyGenericHydro__

#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Physics/GenericHydro.hh"
#include "PyAbstractPhysics.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// PyGenericHydro
//------------------------------------------------------------------------------
template<typename Dimension>
class PyGenericHydro: public PyAbstractPhysics<Dimension, GenericHydro<Dimension>> {
public:
  typedef GenericHydro<Dimension> PhysicsBase;
  typedef PyAbstractPhysics<Dimension, PhysicsBase> PyAP;
  using PyAP::PyAP;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
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
};

}
}

#endif
