//------------------------------------------------------------------------------
// safeInv
//
// Take the inverse of a number, turning it to zero if the argument is 0.0.
//------------------------------------------------------------------------------
#ifndef __Spheral_safevInv__
#define __Spheral_safevInv__

#include "SpheralFunctions.hh"

namespace Spheral {

template<typename Value>
RAJA_HOST_DEVICE
inline
Value
safeInv(const Value& x,
	const double fuzz = 1.0e-30) {
  return x/(x*x + fuzz);
}

template<typename Value>
RAJA_HOST_DEVICE
inline
Value
safeInvVar(const Value& x,
           const double fuzz = 1.0e-30) {
  return sgn(x)/SPHERAL_MAX(fuzz, SPHERAL_ABS(x));
}

}

#endif
