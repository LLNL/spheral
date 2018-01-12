#ifndef __PBGWRAPS_PolyClipperTypes__
#define __PBGWRAPS_PolyClipperTypes__

#include <vector>
#include "Geometry/polyclipper.hh"
#include "PBGWraps/referenceAsPointer.hh"

namespace PolyClipper {
typedef PolyClipper::Plane2d PolyClipperPlane2d;
typedef PolyClipper::Plane3d PolyClipperPlane3d;
}

typedef std::vector<PolyClipper::Plane2d> vector_of_PolyClipperPlane2d;
typedef std::vector<PolyClipper::Plane3d> vector_of_PolyClipperPlane3d;

#endif

