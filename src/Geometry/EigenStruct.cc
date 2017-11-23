//---------------------------------Spheral++----------------------------------//
// EigenStruct, a data structure to hold eigen values and their associated
// eigen vectors.
//
// Created by J. Michael Owen, Thu Dec 28 16:10:05 PST 2000
//----------------------------------------------------------------------------//
#include "EigenStruct.hh"
#include "GeomVector.hh"
#include "GeomTensor.hh"

// Explicit instantiation.
namespace Spheral {
  template struct EigenStruct<1>;
  template struct EigenStruct<2>;
  template struct EigenStruct<3>;
}
