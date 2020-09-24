//---------------------------------Spheral++----------------------------------//
// GeomVector -- Geometric Vector Class.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_hh__
#define __Spheral_GeomVector_hh__

#if defined GEOMMEM_EIGEN
#include "Geometry/GeomVector_eigen.hh"
#elif defined GEOMMEM_ARRAY
#include "Geometry/GeomVector_array.hh"
#else
#include "Geometry/GeomVector_default.hh"
#endif

namespace Spheral {
  template<int nDim> const unsigned GeomVector<nDim>::nDimensions = nDim;
  template<int nDim> const unsigned GeomVector<nDim>::numElements = nDim;
  template<int nDim> const GeomVector<nDim> GeomVector<nDim>::zero = GeomVector<nDim>(0);
  template<int nDim> const GeomVector<nDim> GeomVector<nDim>::one = GeomVector<nDim>(nDim);
}
#endif
