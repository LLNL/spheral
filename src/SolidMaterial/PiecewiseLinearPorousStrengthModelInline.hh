namespace Spheral {

template<typename Dimension>
inline 
std::vector<typename Dimension::Scalar>& 
PiecewiseLinearPorousStrengthModel<Dimension>::
porosityAbscissa(){
  return mPorosityAbscissa;
}

template<typename Dimension>
inline
std::vector<typename Dimension::Scalar>&
PiecewiseLinearPorousStrengthModel<Dimension>::
shearModulusRatios(){
  return mShearModulusRatios;
}

template<typename Dimension>
inline
std::vector<typename Dimension::Scalar>&
PiecewiseLinearPorousStrengthModel<Dimension>::
yieldStrengthRatios(){
  return mYieldStrengthRatios;
}

}



