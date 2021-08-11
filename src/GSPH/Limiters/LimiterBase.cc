#include "LimiterBase.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
LimiterBase<Dimension>::
LimiterBase(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
LimiterBase<Dimension>::
~LimiterBase(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
LimiterBase<Dimension>::
slopeLimiter(const typename Dimension::Scalar x) const {
  return ( x > 0.0 ? 2.0/(1.0+x) * this->fluxLimiter(x) : 0.0 );
}


}