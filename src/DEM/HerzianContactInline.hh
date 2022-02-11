namespace Spheral {

template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
YoungsModulus() const {
  return mYoungsModulus;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
YoungsModulus(typename Dimension::Scalar x) {
  mYoungsModulus = x;
}



template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
restitutionCoefficient() const {
  return mRestitutionCoefficient;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
restitutionCoefficient(typename Dimension::Scalar x) {
  mRestitutionCoefficient = x;
}



template<typename Dimension>
inline
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
beta() const {
  return mBeta;
}
template<typename Dimension>
inline
void
DampedLinearSpring<Dimension>::
beta(typename Dimension::Scalar x) {
  mBeta = x;
}

}