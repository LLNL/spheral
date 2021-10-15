#include "MinModLimiter.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
MinModLimiter<Dimension>::
MinModLimiter():
  LimiterBase<Dimension>(true,true){
}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
MinModLimiter<Dimension>::
~MinModLimiter(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
MinModLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return std::min(1.0, std::max(x, 0.0));
}


}