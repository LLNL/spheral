//------------------------------------------------------------------------------
// This is kinda silly, but we provide an overridable interface to some functors
// here to assist passing functions from Python to C++.
//------------------------------------------------------------------------------
#ifndef __Spheral_Overridable_Functors__
#define __Spheral_Overridable_Functors__

#include "Geometry/Dimension.hh"

namespace Spheral {
namespace PythonBoundFunctors {

// retT F(argT) 
template<typename argT, typename retT>
class SpheralFunctor {
public:
  SpheralFunctor() {};
  virtual ~SpheralFunctor() {};
  virtual retT operator()(const argT x) const { return __call__(x); }
  virtual retT __call__(const argT x) const = 0;
};

}
}

#endif
