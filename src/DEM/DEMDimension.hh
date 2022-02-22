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
      Dim<1>::Scalar cross(const Dim<1>::Scalar v1, const Dim<1>::Scalar v2) const;
  };
  template<>
  class DEMDimension<Dim<2>>{
    public:
      typedef Dim<2>::Scalar AngularVector;
      static const Dim<2>::Scalar zero;

      AngularVector  cross(const Dim<2>::Vector v1, const Dim<2>::Vector v2) const;
      Dim<2>::Vector cross(const AngularVector  v1, const Dim<2>::Vector v2) const;
      Dim<2>::Vector cross(const Dim<2>::Vector v1, const AngularVector  v2) const;
  };
  template<>
  class DEMDimension<Dim<3>>{
    public:
      typedef Dim<3>::Vector AngularVector;
      static const Dim<3>::Vector zero;
      Dim<3>::Vector  cross(const Dim<3>::Vector v1, const Dim<3>::Vector v2) const;
  };
}

#include "DEMDimensionInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> struct DEMDimension;
}

#endif