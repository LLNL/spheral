//------------------------------------------------------------------------------
// Trampoline classes for Integrator interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyIntegrator__
#define __Spheral_PyIntegrator__

#include "Geometry/GeomPlane.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {
namespace IntegratorSpace {

//------------------------------------------------------------------------------
// PyIntegrator
//------------------------------------------------------------------------------
template<typename Dimension, class Base>
class PyIntegrator: public Base {
public:
  using Base::Base;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  virtual void step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      Base,         // Parent class
                      step,         // name of method
                      maxTime,      // arguments
                      state,
                      derivs);
  }

};

}
}

#endif
