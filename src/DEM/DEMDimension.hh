//---------------------------------Spheral++----------------------------------//
// DEM type alias for the particle rotation
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMDimension_hh__
#define __Spheral_DEMDimension_hh__

#include "Geometry/Dimension.hh"

namespace Spheral {
  template<typename Dimension>
  class DEMDimension{};

  template<>
  class DEMDimension<Dim<1>>{
    public:
      typedef Dim<1>::Scalar AngularVector;
      static const Dim<1>::Scalar zero;
  };
  template<>
  class DEMDimension<Dim<2>>{
    public:
      typedef Dim<2>::Scalar AngularVector;
      static const Dim<2>::Scalar zero;
  };
  template<>
  class DEMDimension<Dim<3>>{
    public:
      typedef Dim<3>::Vector AngularVector;
      static const Dim<3>::Vector zero;
  };
}
#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> struct DEMDimension;
}

#endif