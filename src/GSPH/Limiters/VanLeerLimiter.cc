//---------------------------------Spheral++----------------------------------//
// VanLeerLimiter
//   Van Leer, B. (1979), "Towards the ultimate conservative difference scheme 
//   V. A second order sequel to Godunov's method", J. Comput. Phys., 32 (1): 
//   101–136
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "VanLeerLimiter.hh"

namespace Spheral {

//----------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------
template<typename Dimension>
VanLeerLimiter<Dimension>::
VanLeerLimiter():
  LimiterBase<Dimension>(true,true){
}

//----------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------
template<typename Dimension>
VanLeerLimiter<Dimension>::
~VanLeerLimiter(){}

//----------------------------------------------------------------------------
// slope limiter
//----------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
VanLeerLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return 2.0*x/(1.0 + x) ;
}


}