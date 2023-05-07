namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
CircularFinitePlane<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
CircularFinitePlane<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
CircularFinitePlane<Dimension>::
normal() const {
  return mNormal;
}

template<typename Dimension>
inline
void
CircularFinitePlane<Dimension>::
normal(const typename Dimension::Vector& value) {
  mNormal = value;
}

template<typename Dimension>
inline
const typename Dimension::Scalar&
CircularFinitePlane<Dimension>::
extent() const {
  return mExtent;
}

template<typename Dimension>
inline
void
CircularFinitePlane<Dimension>::
extent(const typename Dimension::Scalar& value) {
  mExtent=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
CircularFinitePlane<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
CircularFinitePlane<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}