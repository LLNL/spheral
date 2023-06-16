namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
center() const {
  return mCenter;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
center(const typename Dimension::Vector& value) {
  mCenter=value;
}

template<typename Dimension>
inline
typename Dimension::Scalar
SphereSolidBoundary<Dimension>::
radius() const {
  return mRadius;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
radius(typename Dimension::Scalar value) {
  mRadius=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
clipPoint() const {
  return mClipPoint;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
clipPoint(const typename Dimension::Vector& value) {
  mClipPoint=value;
  this->setClipIntersectionRadius();
}


template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
clipAxis() const {
  return mClipAxis;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
clipAxis(const typename Dimension::Vector& value) {
  mClipAxis=value.unitVector();
  this->setClipIntersectionRadius();
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}