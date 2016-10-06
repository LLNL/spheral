//---------------------------------Spheral++----------------------------------//
// hashes
//
// A collection of simple useful hashes.
//
// Created by JMO, Fri Jan  7 21:56:22 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_hashes__
#define __Spheral_hashes__

namespace Spheral {

//------------------------------------------------------------------------------
// Combine two unsigneds to a single unsigned long long.
//------------------------------------------------------------------------------
inline
uint64_t
hashUnsignedPair(uint32_t a, uint32_t b) {
  if (a > b) std::swap(a, b);
  uint64_t result = (uint64_t) a << 32 | b;
  return result;
}

}

#endif

