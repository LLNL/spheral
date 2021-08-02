#include "SlopeLimiterBase.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
SlopeLimiterBase<Dimension>::
SlopeLimiterBase(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
SlopeLimiterBase<Dimension>::
~SlopeLimiterBase(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
SlopeLimiterBase<Dimension>::
slopeLimiter(const typename Dimension::Scalar x) const {
  return ( x > 0.0 ? 2.0/(1.0+x) * this->fluxLimiter(x) : 0.0 );
}


}