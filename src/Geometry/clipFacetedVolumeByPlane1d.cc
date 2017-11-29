//---------------------------------Spheral++----------------------------------//
// Clip a faceted volume (box, polygon, or polyhedron) by a plane in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) >= 0
// is retained.
//
// The algorithms herein are roughly based on the approach outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
// though I think not exactly the same.
//
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//

#include "clipFacetedVolumeByPlane.hh"

namespace Spheral {

void clipFacetedVolumeByPlane(Box1d& poly, const GeomPlane<Dim<1>>& plane) {

  // Useful types.
  typedef Dim<1>::FacetedVolume FacetedVolume;
  typedef Dim<1>::Vector Vector;

  const auto boxmin = poly.xmin();
  const auto boxmax = poly.xmax();
  auto center = poly.center();
  auto extent = poly.extent();

  if (plane.compare(boxmin) >= 0 and plane.compare(boxmax) >= 0) {
    // The input box is entirely above the plane, so is unclipped.
    return;

  } else if (plane.compare(boxmin) < 0 and plane.compare(boxmax) < 0) {
    // The input box is entirely below the plane, and therefore is entirely rejected.
    center = Vector::zero;
    extent = 0.0;

  } else {
    // The plane splits our box somewhere, so we actually need to clip it.
    const auto& p = plane.point();
    const auto& phat = plane.normal();
    CHECK(poly.contains(p));
    CHECK(phat.x() == 1.0 or phat.x() == -1.0);
    if (phat.x() > 0.0) {
      center = 0.5*(boxmin + p);
      extent = center.x() - boxmin.x();
    } else {
      center = 0.5*(boxmax + p);
      extent = boxmax.x() - center.x();
    }
    CHECK(extent >= 0.0);

    poly.center(center);
    poly.extent(extent);
  }
}

}
