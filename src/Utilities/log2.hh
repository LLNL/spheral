//---------------------------------Spheral++----------------------------------//
// log2 -- various methods for computing base 2 logarithms.
//
// Created by JMO, Tue Feb  9 21:44:34 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_log2__
#define __Spheral_log2__

#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Find the integer log2 of an integer (unsigned of course).
// This is a dumb direct method, not fast!
//------------------------------------------------------------------------------
inline
unsigned
log2(unsigned x) {
  REQUIRE(x > 0);
  unsigned result = 1U;
  while (result < x) result <<= 1;
  result >>= 1;
  ENSURE((1 << result) <= x and x <  (1 << (result + 1)));
  return result;
}

//------------------------------------------------------------------------------
// Find the integer log2 of an integer (unsigned of course).
// This is based on code found at 
// http://graphics.stanford.edu/~seander/bithacks.html
//------------------------------------------------------------------------------
inline
unsigned
log2(uint32_t x) {
  REQUIRE(x > 0);
  static const int MultiplyDeBruijnBitPosition[32] = {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
  };

  x |= x >> 1; // first round down to one less than a power of 2 
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;

  return MultiplyDeBruijnBitPosition[(uint32_t)(x * 0x07C4ACDDU) >> 27];
}

}

#endif

