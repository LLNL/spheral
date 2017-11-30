//---------------------------------Spheral++----------------------------------//
// Clip a faceted volume (box, polygon, or polyhedron) by a plane in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) <= 0
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
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;

namespace {

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
segmentPlaneIntersection(const Dim<2>::Vector& a,       // line-segment begin
                         const Dim<2>::Vector& b,       // line-segment end
                         const Dim<2>::Vector& p,       // point in plane
                         const Dim<2>::Vector& phat) {  // plane unit normal

  const auto ab = b - a;
  const auto abhat = ab.unitVector();
  CHECK2(std::abs(abhat.dot(phat)) > 0.0, (abhat.dot(phat)) << " " << a << " " << b << " " << abhat << " " << phat);
  const auto s = (p - a).dot(phat)/(abhat.dot(phat));
  CHECK(s >= 0.0 and s <= ab.magnitude());
  const auto result = a + s*abhat;
  CHECK2(fuzzyEqual((result - p).dot(phat), 0.0, 1.0e-10),
         a << " " << b << " " << s << " " << result << " " << (result - p).dot(phat));
  return result;
}

//------------------------------------------------------------------------------
// Compare a plane and point (our built-in plane one has some issues).
//------------------------------------------------------------------------------
inline
int compare(const Dim<2>::Vector& planePoint,
            const Dim<2>::Vector& planeNormal,
            const Dim<2>::Vector& point) {
  return sgn0(planeNormal.dot(point - planePoint));
}

//------------------------------------------------------------------------------
// Extract the ordered ring of vertices making up a polygon.
//------------------------------------------------------------------------------
vector<unsigned> orderedPolygonRing(const GeomPolygon& poly) {
  
  const auto& facets = poly.facets();
  const auto  nf = facets.size();
  CHECK(nf >= 3);

  // Kickstart with the first facet.
  vector<unsigned> result;
  result.push_back(facets[0].ipoint1());
  result.push_back(facets[0].ipoint2());

  // Walk the facets, inserting in order.
  for (unsigned k = 1; k < nf; ++k) {
    if (result.back() == facets[k].ipoint1()) {
      // The most common case where the facets themselves have been ordered.
      result.push_back(facets[k].ipoint2());
    } else {
      // See if our starting point is already in the result.
      auto itr = find(result.begin(), result.end(), facets[k].ipoint1());
      if (itr != result.end()) {
        // Our starting point is already there, so insert our final point in order.
        if ( (itr+1) == result.end() or
            *(itr+1) != facets[k].ipoint2()) {
          result.insert(itr, facets[k].ipoint2());
        }
      } else {
        // The points of this facet are unknown, so huck em in.
        result.push_back(facets[k].ipoint1());
        result.push_back(facets[k].ipoint2());
      }
    }
  }

  // The first and last points should now be redundant, so clip that last one.
  CHECK(result.size() >= 4);
  CHECK(result.front() == result.back());
  result.pop_back();

  // // BLAGO
  // const auto& vertices = poly.vertices();
  // cerr << "Unique ring: " << endl
  //      << " ---> ";
  // for (auto k = 0; k < result.size(); ++k) cerr << result[k] << " ";
  // cerr << endl
  //      << " ---> ";
  // for (auto k = 0; k < result.size(); ++k) cerr << vertices[result[k]] << " ";
  // cerr << endl;
  // // BLAGO

  // That's it.
  return result;
}

}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
void clipFacetedVolumeByPlanes(GeomPolygon& poly, 
                               const std::vector<GeomPlane<Dim<2>>>& planes) {

  // Useful types.
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::Vector Vector;

  // The following logic depends on the vertices of the polygon already being arranged 
  // in counter-clockwise order.
  auto ring0 = poly.vertices();    // Note this is a copy!
  auto nring0 = ring0.size();

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not ring0.empty()) {
    const auto& plane = planes[kplane++];

    // Check if the polygon is entirely clear of the plane (above or below).
    auto above = true;
    auto below = true;
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();
    {
      auto k = 0;
      while (k < nring0 and (above or below)) {
        const auto vcomp = compare(p0, phat, ring0[k]);
        if (vcomp >= 0) {
          below = false;
        } else if (vcomp == -1) {
          above = false;
        }
        ++k;
      }
    }
    CHECK(not (above and below));

    // Did we get a simple case?
    if (below) {
      // Polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      ring0.clear();

    } else if (not above) {

      // The plane passes somewhere through the polygon, so we need to walk the ring and modify it.
      vector<Vector> ring1;
      for (auto k = 0; k < nring0; ++k) {
        const auto& v0 = ring0[k];
        const auto& v1 = ring0[(k + 1) % nring0];
        const auto comp0 = compare(p0, phat, v0);
        const auto comp1 = compare(p0, phat, v1);
        if (comp0 >= 0) {
          ring1.push_back(v0);
          if (comp1 < 0) {
            const auto pintersect = segmentPlaneIntersection(v0, v1, p0, phat);
            ring1.push_back(pintersect);
          }
        } else if (comp1 > 0) {
          const auto pintersect = segmentPlaneIntersection(v0, v1, p0, phat);
          ring1.push_back(pintersect);
        }
      }
      ring0 = ring1;
      nring0 = ring0.size();
    }
  }

  // Now rebuild the polygon, and we're done.
  if (ring0.empty()) {
    poly = GeomPolygon();
  } else {
    CHECK(nring0 >= 3);
    vector<vector<unsigned>> facets(nring0, vector<unsigned>(2));
    for (auto k = 0; k < nring0; ++k) {
      facets[k][0] = k;
      facets[k][1] = (k + 1) % nring0;
    }
    poly = GeomPolygon(ring0, facets);
  }
}

}
