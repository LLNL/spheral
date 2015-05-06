//---------------------------------Spheral++----------------------------------//
// GeomVector_fwd -- forward declaration for Geometric Vector Class.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_fwd_hh__
#define __Spheral_GeomVector_fwd_hh__

namespace Spheral {
#ifdef GEOMMEM
  template<int nDim, bool ownMemory = true> class GeomVector;
#else
  template<int nDim> class GeomVector;
#endif
}

#endif
