namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularFinitePlane<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
RectangularFinitePlane<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularFinitePlane<Dimension>::
extent() const {
  return mExtent;
}

template<typename Dimension>
inline
void
RectangularFinitePlane<Dimension>::
extent(const typename Dimension::Vector& value) {
  mExtent=value;
}

template<typename Dimension>
inline
const typename Dimension::Tensor&
RectangularFinitePlane<Dimension>::
basis() const {
  return mBasis;
}

template<typename Dimension>
inline
void
RectangularFinitePlane<Dimension>::
basis(const typename Dimension::Tensor& value) {
  mBasis=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularFinitePlane<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
RectangularFinitePlane<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}