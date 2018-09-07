namespace Spheral {

//------------------------------------------------------------------------------
// Access the entrance plane.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>&
PlanarBoundary<Dimension>::enterPlane() const {
  return mEnterPlane;
}

template<typename Dimension>
inline
void
PlanarBoundary<Dimension>::
setEnterPlane(const GeomPlane<Dimension>& enterPlane) {
  mEnterPlane = enterPlane;
}

//------------------------------------------------------------------------------
// Access the exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>&
PlanarBoundary<Dimension>::exitPlane() const {
  return mExitPlane;
}

template<typename Dimension>
inline
void
PlanarBoundary<Dimension>::
setExitPlane(const GeomPlane<Dimension>& exitPlane) {
  mExitPlane = exitPlane;
}

}
