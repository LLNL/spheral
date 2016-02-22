//---------------------------------Spheral++----------------------------------//
// planarReflectingOperator
// Generate the reflection operator for the given plane.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_planarReflectingOperator__
#define __Spheral_planarReflectingOperator__

namespace Spheral {

template<typename Plane>
inline
typename Plane::Tensor
planarReflectingOperator(const Plane& plane) {
  typedef typename Plane::Vector Vector;
  typedef typename Plane::Tensor Tensor;
  return Tensor::one - 2.0*plane.normal().selfdyad();
}

}

#endif
