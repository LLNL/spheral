namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
PlanarWall<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
PlanarWall<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
PlanarWall<Dimension>::
normal() const {
  return mNormal;
}

template<typename Dimension>
inline
void
PlanarWall<Dimension>::
normal(const typename Dimension::Vector& value)  {
  mNormal=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
PlanarWall<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
PlanarWall<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}