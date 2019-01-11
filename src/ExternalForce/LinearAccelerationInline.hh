#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the acceleration parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearAcceleration<Dimension>::
a0() const {
  return ma0;
}

template<typename Dimension>
inline
typename Dimension::Scalar
LinearAcceleration<Dimension>::
aslope() const {
  return maslope;
}

}
