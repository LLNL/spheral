//---------------------------------Spheral++----------------------------------//
// DEM type alias for the particle rotation
//----------------------------------------------------------------------------//
namespace Spheral{
    
  // cross product
  // -------------------------------------------------------------------------/
  inline
  typename DEMDimension<Dim<1>>::AngularVector  
  DEMDimension<Dim<1>>::
  cross(const Dim<1>::Vector v1, const Dim<1>::Vector v2){
    return 0.0;
  }

  inline
  typename Dim<1>::Vector  
  DEMDimension<Dim<1>>::
  cross(const DEMDimension<Dim<1>>::AngularVector  v1, const Dim<1>::Vector v2){
    return Dim<1>::Vector(0.0);
  }

  inline
  typename Dim<1>::Vector  
  DEMDimension<Dim<1>>::
  cross(const Dim<1>::Vector v1, const DEMDimension<Dim<1>>::AngularVector  v2){
    return Dim<1>::Vector(0.0);
  }

  inline
  typename DEMDimension<Dim<2>>::AngularVector  
  DEMDimension<Dim<2>>::
  cross(const Dim<2>::Vector v1, const Dim<2>::Vector v2){
    return 0.0;
  }

  inline
  typename Dim<2>::Vector  
  DEMDimension<Dim<2>>::
  cross(const DEMDimension<Dim<2>>::AngularVector  v1, const Dim<2>::Vector v2){
    return Dim<2>::Vector(-v1*v2[1], v1*v2[0]);
  }

  inline
  typename Dim<2>::Vector  
  DEMDimension<Dim<2>>::
  cross(const Dim<2>::Vector v1, const DEMDimension<Dim<2>>::AngularVector  v2){
    return Dim<2>::Vector(v1[0]*v2, -v1[1]*v2);
  }

  inline
  typename Dim<3>::Vector  
  DEMDimension<Dim<3>>::
  cross(const Dim<3>::Vector v1, const Dim<3>::Vector v2){
    return Dim<3>::Vector(v1[0]*v2[2] - v1[2]*v2[1],
                          v1[2]*v2[0] - v1[0]*v2[2],
                          v1[0]*v2[1] - v1[1]*v2[0]);
  }
  
  // dot product this is only nonzero for 3d
  // -------------------------------------------------------------------------/
  inline
  typename Dim<1>::Scalar
  DEMDimension<Dim<1>>::
  dot(const Dim<1>::Vector, const DEMDimension<Dim<1>>::AngularVector){
    return 0.0;
  }
  inline
  typename Dim<1>::Scalar
  DEMDimension<Dim<1>>::
  dot(const DEMDimension<Dim<1>>::AngularVector,const Dim<1>::Vector){
    return 0.0;
  }
  inline
  typename Dim<2>::Scalar
  DEMDimension<Dim<2>>::
  dot(const Dim<2>::Vector, const DEMDimension<Dim<2>>::AngularVector){
    return 0.0;
  }
  inline
  typename Dim<2>::Scalar
  DEMDimension<Dim<2>>::
  dot(const DEMDimension<Dim<2>>::AngularVector,const Dim<2>::Vector){
    return 0.0;
  }
  inline
  typename Dim<3>::Scalar
  DEMDimension<Dim<3>>::
  dot(const Dim<3>::Vector v1, const Dim<3>::Vector v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }
}
