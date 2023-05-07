namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlane<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
InfinitePlane<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlane<Dimension>::
normal() const {
  return mNormal;
}

template<typename Dimension>
inline
void
InfinitePlane<Dimension>::
normal(const typename Dimension::Vector& value) {
  mNormal = value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlane<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
InfinitePlane<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}