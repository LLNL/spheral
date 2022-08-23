namespace Spheral {


inline
const std::vector<Dim<3>::Vector>& 
ApproximatePolyhedralGravityModel::values() const {
  return mValues;
}  

inline
const std::vector<Dim<3>::Vector>& 
ApproximatePolyhedralGravityModel::quadraturePoints() const {
  return mQuadraturePoints;
}   

inline
const std::vector<Dim<3>::Scalar>& 
ApproximatePolyhedralGravityModel::resolutions() const {
  return mResolutions;
}   

inline
Dim<3>::Scalar
ApproximatePolyhedralGravityModel::numQuadraturePoints() const {
  return mNumQuadraturePoints;
}   


}
