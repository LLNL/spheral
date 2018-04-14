//------------------------------------------------------------------------------
// Trampoline classes for Integrator abstract interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyAbstractIntegrator__
#define __Spheral_PyAbstractIntegrator__

#include "Geometry/GeomPlane.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {
namespace IntegratorSpace {

//------------------------------------------------------------------------------
// PyAbstractIntegrator
//------------------------------------------------------------------------------
template<typename Dimension, typename Base>
class PyAbstractIntegrator: public Base {
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
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           Base,         // Parent class
                           step,         // name of method
                           maxTime,      // arguments
                           state,
                           derivs
      );
  }

    virtual void step(Scalar maxTime) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      Base,         // Parent class
                      step,         // name of method
                      maxTime);     // arguments
  }

  virtual Scalar selectDt(const Scalar dtMin, 
                          const Scalar dtMax,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs) const override {
    PYBIND11_OVERLOAD(Scalar,       // Return type
                      Base,         // Parent class
                      selectDt,     // name of method
                      dtMin,        // arguments
                      dtMax,
                      state,
                      derivs);
  }

  virtual void preStepInitialize(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,              // Return type
                      Base,              // Parent class
                      preStepInitialize, // name of method
                      state, // arguments
                      derivs);
  }

  virtual void initializeDerivatives(const double t,
                                     const double dt,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      Base,         // Parent class
                      initializeDerivatives,         // name of method
                      t,            // arguments
                      dt,
                      state,
                      derivs);
  }

  virtual void postStepFinalize(const double t,
                                const double dt,
                                State<Dimension>& state,
                                StateDerivatives<Dimension>& derivs) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      Base,         // Parent class
                      postStepFinalize,         // name of method
                      t,            // arguments
                      dt,
                      state,
                      derivs);
  }

  virtual void advance(Scalar goalTime) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      Base,         // Parent class
                      advance,      // name of method
                      goalTime);    // arguments
  }

};

}
}

#endif
