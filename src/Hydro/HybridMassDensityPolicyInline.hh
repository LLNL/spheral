namespace Spheral {

//------------------------------------------------------------------------------
// Weight scale value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
HybridMassDensityPolicy<Dimension>::
weightScale() const {
  return sqrt(mWeightScale2);
}

//------------------------------------------------------------------------------
// Min value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
HybridMassDensityPolicy<Dimension>::
minValue() const {
  return mMinValue;
}

//------------------------------------------------------------------------------
// Max value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
HybridMassDensityPolicy<Dimension>::
maxValue() const {
  return mMaxValue;
}

}
