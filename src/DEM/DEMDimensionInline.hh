//---------------------------------Spheral++----------------------------------//
// DEM type alias for the particle rotation
//----------------------------------------------------------------------------//
namespace Spheral{
    
  inline
  typename Dim<1>::Scalar
  DEMDimension<Dim<1>>::
  cross(const Dim<1>::Scalar v1, const Dim<1>::Scalar v2) const {
    return 0;
  }

  inline
  typename DEMDimension<Dim<2>>::AngularVector  
  DEMDimension<Dim<2>>::
  cross(const Dim<2>::Vector v1, const Dim<2>::Vector v2) const {
    return (v1[0]*v2[1] - v1[1]*v2[0]);
  }

  inline
  typename Dim<2>::Vector  
  DEMDimension<Dim<2>>::
  cross(const DEMDimension<Dim<2>>::AngularVector  v1, const Dim<2>::Vector v2) const {
    return Dim<2>::Vector(v1*v2[1], - v1*v2[0]);
  }

  inline
  typename Dim<2>::Vector  
  DEMDimension<Dim<2>>::
  cross(const Dim<2>::Vector v1, const DEMDimension<Dim<2>>::AngularVector  v2) const {
    return Dim<2>::Vector(-v1[0]*v2, v1[1]*v2);
  }

  inline
  typename Dim<3>::Vector  
  DEMDimension<Dim<3>>::
  cross(const Dim<3>::Vector v1, const Dim<3>::Vector v2) const {
    return Dim<3>::Vector(v1[0]*v2[2] - v1[2]*v2[1],
                          v1[2]*v2[0] - v1[0]*v2[2],
                          v1[0]*v2[1] - v1[1]*v2[0]);
  }
  
}
