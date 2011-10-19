#include "cdebug.hh"

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Access the entrance plane.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>&
PlanarBoundary<Dimension>::enterPlane() const {
  cdebug << "PlanarBoundary::enterPlane(): " << this << std::endl;
  return mEnterPlane;
}

template<typename Dimension>
inline
void
PlanarBoundary<Dimension>::
setEnterPlane(const GeomPlane<Dimension>& enterPlane) {
  cdebug << "PlanarBoundary::setEnterPlane(enterPlane): " << this << std::endl;
  mEnterPlane = enterPlane;
}

//------------------------------------------------------------------------------
// Access the exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>&
PlanarBoundary<Dimension>::exitPlane() const {
  cdebug << "PlanarBoundary::exitPlane(): " << this << std::endl;
  return mExitPlane;
}

template<typename Dimension>
inline
void
PlanarBoundary<Dimension>::
setExitPlane(const GeomPlane<Dimension>& exitPlane) {
  cdebug << "PlanarBoundary::setExitPlane(exitPlane): " << this << std::endl;
  mExitPlane = exitPlane;
}

}
}
