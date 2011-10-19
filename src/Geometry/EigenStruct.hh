//---------------------------------Spheral++----------------------------------//
// EigenStruct, a data structure to hold eigen values and their associated
// eigen vectors.
//
// Created by J. Michael Owen, Thu Dec 28 16:10:05 PST 2000
//----------------------------------------------------------------------------//
#ifndef EigenStruct_HH
#define EigenStruct_HH

#include <iostream>

namespace Spheral {

template<int nDim> class GeomVector;
template<int nDim> class GeomTensor;

// EigenStruct
template<int nDim>
struct EigenStruct {
  GeomVector<nDim> eigenValues;
  GeomTensor<nDim> eigenVectors;

  EigenStruct() {}
  EigenStruct(const EigenStruct& rhs):
    eigenValues(rhs.eigenValues),
    eigenVectors(rhs.eigenVectors) {
  }
  EigenStruct& operator=(const EigenStruct& rhs) {
    if (this != &rhs) {
      eigenValues = rhs.eigenValues;
      eigenVectors = rhs.eigenVectors;
    }
    return *this;
  }
};

// Forward declarations.
template<int nDim> std::istream& operator>>(::std::istream is, Spheral::EigenStruct<nDim>& eigen);
template<int nDim> std::ostream& operator<<(::std::ostream os, Spheral::EigenStruct<nDim>& eigen);

}

// Include these here to avoid compiler complaints incomplete types.
#include "GeomVector.hh"
#include "GeomTensor.hh"

#ifndef __GCCXML__
#include "EigenStructInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<int nDim> struct EigenStruct;
}

#endif
