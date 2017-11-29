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

#include "clipFacetedVolumeByPlanes.hh"
#include <limits>

namespace Spheral {

void clipFacetedVolumeByPlanes(Box1d& poly, 
                               const std::vector<GeomPlane<Dim<1>>>& planes) {

  // Useful types.
  typedef Dim<1>::FacetedVolume FacetedVolume;
  typedef Dim<1>::Vector Vector;

  // Look for the min and max allowed range across the planes.
  auto clipmin = std::numeric_limits<double>::max();
  auto clipmax = std::numeric_limits<double>::lowest();
  for (const auto& plane: planes) {
    if (plane.normal().x() > 0.0) {
      clipmin = std::max(clipmin, plane.point().x());
    } else {
      clipmax = std::min(clipmax, plane.point().x());
    }
  }

  auto boxmin = poly.xmin().x();
  auto boxmax = poly.xmax().x();
  auto center = poly.center().x();
  auto extent = poly.extent();

  if (clipmin >= clipmax or
      clipmin >= boxmax or
      clipmax <= boxmin) {
    // The box is entirely excluded, so empty.
    center = 0.0;
    extent = 0.0;

  } else {

    boxmin = std::max(boxmin, clipmin);
    boxmax = std::min(boxmax, clipmax);
    CHECK(boxmax - boxmin >= 0.0);
    center = 0.5*(boxmin + boxmax);
    extent = 0.5*(boxmax - boxmin);

  }

  // Update the box.
  poly.center(Vector(center));
  poly.extent(extent);
}

}
