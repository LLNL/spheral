//---------------------------------Spheral++----------------------------------//
// OspreLimiter 
//   Waterson, N.P.; Deconinck, H. (1995), "A unified approach to the design 
//   and application of bounded higher-order convection schemes," Journal
//   of Computational Physics, 224 (1): 182-207. 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "OspreLimiter.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
OspreLimiter<Dimension>::
OspreLimiter():
  LimiterBase<Dimension>(true,true){
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
OspreLimiter<Dimension>::
~OspreLimiter(){}

//------------------------------------------------------------------------------
// slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
OspreLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return 1.5 * (x*x + x) / (x*x + x + 1.0) ;
}


}