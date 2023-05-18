namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
FiniteCylinder<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
FiniteCylinder<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
FiniteCylinder<Dimension>::
axis() const {
  return mAxis;
}

template<typename Dimension>
inline
void
FiniteCylinder<Dimension>::
axis(const typename Dimension::Vector& value) {
  mAxis=value;
}

template<typename Dimension>
inline
typename Dimension::Scalar
FiniteCylinder<Dimension>::
length() const {
  return mLength;
}

template<typename Dimension>
inline
void
FiniteCylinder<Dimension>::
length(typename Dimension::Scalar value) {
  mLength=value;
}

template<typename Dimension>
inline
typename Dimension::Scalar
FiniteCylinder<Dimension>::
radius() const {
  return mRadius;
}

template<typename Dimension>
inline
void
FiniteCylinder<Dimension>::
radius(typename Dimension::Scalar value) {
  mRadius=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
FiniteCylinder<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
FiniteCylinder<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}