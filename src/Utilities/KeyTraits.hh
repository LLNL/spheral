//---------------------------------Spheral++----------------------------------//
// KeyTraits
//
// Encapsulate how we think about keys for the space filling curves.
// This is specialized assuming we're working with 64 bit uint64_t.
//
// Created by JMO, Wed Apr  9 13:21:29 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_Utilities_KeyTraits__
#define __Spheral_Utilities_KeyTraits__

#include <stdint.h>

namespace Spheral {

struct KeyTraits {
  typedef uint64_t Key;
  static const uint32_t numbits;
  static const uint32_t numbits1d;
  static const Key zero;
  static const Key one;
  static const Key two;
  static const Key maxKey1d;
  static const Key maxKey;
};

}

#endif
