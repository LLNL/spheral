#include "VanLeerLimiter.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
VanLeerLimiter<Dimension>::
VanLeerLimiter(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
VanLeerLimiter<Dimension>::
~VanLeerLimiter(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
VanLeerLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return 2.0*x/(1.0 + x) ;
}


}