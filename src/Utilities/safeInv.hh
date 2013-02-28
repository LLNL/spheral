//------------------------------------------------------------------------------
// safeInv
//
// Take the inverse of a number, turning it to zero if the argument is 0.0.
//------------------------------------------------------------------------------
#ifndef __Spheral_safevInv__
#define __Spheral_safevInv__

template<typename Value>
inline
Value
safeInv(const Value& x,
	const double fuzz = 1.0e-30) {
  return x/(x*x + fuzz);
}

#endif
