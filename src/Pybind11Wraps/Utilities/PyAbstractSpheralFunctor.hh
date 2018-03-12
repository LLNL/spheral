//------------------------------------------------------------------------------
// Trampoline classes for SpheralFunctor abstract interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyAbstractSpheralFunctor__
#define __Spheral_PyAbstractSpheralFunctor__

#include "Geometry/GeomPlane.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {

//------------------------------------------------------------------------------
// PyAbstractSpheralFunctor
//------------------------------------------------------------------------------
template<typename argT, typename retT, typename Base>
class PyAbstractSpheralFunctor: public Base {
public:
  using Base::Base;  // inherit constructors

  virtual retT __call__(const argT x) const override {
    PYBIND11_OVERLOAD_PURE(retT,         // Return type
                           Base,         // Parent class
                           __call__,     // name of method
                           x);           // arguments
  }
};

}

#endif
