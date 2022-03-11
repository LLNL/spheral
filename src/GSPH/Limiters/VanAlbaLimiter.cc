//---------------------------------Spheral++----------------------------------//
// VanAlbaLimiter 
//   Van Albada, G.D.; Van Leer, B.; Roberts, W.W. (1982), 
//   "A comparative study of computational methods in cosmic gas dynamics", 
//   Astronomy and Astrophysics, 108 (1): 76–84
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "VanAlbaLimiter.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
VanAlbaLimiter<Dimension>::
VanAlbaLimiter():
  LimiterBase<Dimension>(true,true){
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
VanAlbaLimiter<Dimension>::
~VanAlbaLimiter(){}

//------------------------------------------------------------------------------
// slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
VanAlbaLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return (x*x + x) / (x*x + 1.0);
}


}