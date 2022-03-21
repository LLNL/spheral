//---------------------------------Spheral++----------------------------------//
// LimiterBase -- base class for scalar slope limiters
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "LimiterBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
LimiterBase<Dimension>::
LimiterBase(bool TVD,
            bool symmetric):
  mTVD(TVD),
  mSymmetric(symmetric){

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
LimiterBase<Dimension>::
~LimiterBase(){}

//------------------------------------------------------------------------------
// slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LimiterBase<Dimension>::
slopeLimiter(const typename Dimension::Scalar x) const {
  return ( x > 0.0 ? 2.0/(1.0+x) * this->fluxLimiter(x) : 0.0 );
}


//------------------------------------------------------------------------------
// is this TVD?
//------------------------------------------------------------------------------
template<typename Dimension>
bool
LimiterBase<Dimension>::
isTVD() const {
  return mTVD;
}

//------------------------------------------------------------------------------
// is this Symmetric?
//------------------------------------------------------------------------------
template<typename Dimension>
bool
LimiterBase<Dimension>::
isSymmetric() const {
  return mSymmetric;
}

}