namespace Spheral {

template<typename Dimension>
inline
typename Dimension::Scalar
HerzianDEM<Dimension>::
YoungsModulus() const {
  return mYoungsModulus;
}
template<typename Dimension>
inline
void
HerzianDEM<Dimension>::
YoungsModulus(typename Dimension::Scalar x) {
  mYoungsModulus = x;
}



template<typename Dimension>
inline
typename Dimension::Scalar
HerzianDEM<Dimension>::
restitutionCoefficient() const {
  return mRestitutionCoefficient;
}
template<typename Dimension>
inline
void
HerzianDEM<Dimension>::
restitutionCoefficient(typename Dimension::Scalar x) {
  mRestitutionCoefficient = x;
}



template<typename Dimension>
inline
typename Dimension::Scalar
HerzianDEM<Dimension>::
beta() const {
  return mBeta;
}

template<typename Dimension>
inline
void
HerzianDEM<Dimension>::
beta(typename Dimension::Scalar x) {
  mBeta = x;
}

}