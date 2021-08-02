#include "OspreLimiter.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
OspreLimiter<Dimension>::
OspreLimiter(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
OspreLimiter<Dimension>::
~OspreLimiter(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
OspreLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return 1.5 * (x*x + x) / (x*x + x + 1.0) ;
}


}