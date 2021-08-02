#include "VanAlbaLimiter.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
VanAlbaLimiter<Dimension>::
VanAlbaLimiter(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
VanAlbaLimiter<Dimension>::
~VanAlbaLimiter(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
VanAlbaLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return (x*x + x) / (x*x + 1.0);
}


}