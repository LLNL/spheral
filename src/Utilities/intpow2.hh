//---------------------------------Spheral++----------------------------------//
// intpow2
//
// Compute 2^n, where n is an integer in the range [0, nmax], where nmax is
// the number of bits in an unsigned long long.
//
// Created by JMO, Thu Jan  7 13:09:14 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_intpow2__
#define __Spheral_intpow2__

#include <limits>
#include "Utilities/DBC.hh"

namespace Spheral {

inline
unsigned long long
intpow2(const unsigned n) {
  REQUIRE(n <= std::numeric_limits<unsigned long long>::digits);
  return 1ULL << n;
}

}

#endif
