#include "SuperbeeLimiter.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
SuperbeeLimiter<Dimension>::
SuperbeeLimiter():
  LimiterBase<Dimension>(true,true){
}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
SuperbeeLimiter<Dimension>::
~SuperbeeLimiter(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
SuperbeeLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return std::max(std::min(2.0*x, 1.0), std::min(x, 2.0));
}


}