//---------------------------------PolyClipper--------------------------------//
// Clip a faceted volume (polygon or polyhedron) by a set of planes in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) >= 0
// is retained.
//
// The algorithms herein are based on R3D as outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
//
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//

#include "polyclipper.hh"
#include "polyclipper_utilities.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"
#include "Utilities/removeElements.hh"

#include <list>
#include <map>
#include <iostream>
#include <iterator>
#include <algorithm>

// Declare the timers.
extern Timer TIME_PC2d_convertto;
extern Timer TIME_PC2d_convertfrom;
extern Timer TIME_PC2d_moments;
extern Timer TIME_PC2d_clip;
extern Timer TIME_PC2d_planes;
extern Timer TIME_PC2d_checkverts;
extern Timer TIME_PC2d_insertverts;
extern Timer TIME_PC2d_hanging;
extern Timer TIME_PC2d_compress;

namespace PolyClipper {

using namespace std;

namespace {    // anonymous methods

//------------------------------------------------------------------------------
// Compare a plane and point.
//------------------------------------------------------------------------------
inline
int compare(const Plane2d& plane,
            const Spheral::Dim<2>::Vector& point) {
  const auto sgndist = plane.dist + plane.normal.dot(point);
  if (std::abs(sgndist) < 1.0e-10) return 0;
  return sgn0(sgndist);
}

//------------------------------------------------------------------------------
// Compare a plane and a box (defined by it's min/max coordinates).
//   -1 ==> box below plane
//    0 ==> plane cuts through box
//    1 ==> box above plane
//------------------------------------------------------------------------------
inline
int compare(const Plane2d& plane,
            const double xmin,
            const double ymin,
            const double xmax,
            const double ymax) {
  typedef Spheral::Dim<2>::Vector Vector;
  const auto c1 = compare(plane, Vector(xmin, ymin));
  const auto c2 = compare(plane, Vector(xmax, ymin));
  const auto c3 = compare(plane, Vector(xmax, ymax));
  const auto c4 = compare(plane, Vector(xmin, ymax));
  const auto cmin = std::min(c1, std::min(c2, std::min(c3, c4)));
  const auto cmax = std::max(c1, std::max(c2, std::max(c3, c4)));
  if (cmin >= 0) {
    return  1;
  } else if (cmax <= 0) {
    return -1;
  } else {
    return  0;
  }
}

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
inline
Spheral::Dim<2>::Vector
segmentPlaneIntersection(const Spheral::Dim<2>::Vector& a,       // line-segment begin
                         const Spheral::Dim<2>::Vector& b,       // line-segment end
                         const Plane2d& plane) {                 // plane
  const auto asgndist = plane.dist + plane.normal.dot(a);
  const auto bsgndist = plane.dist + plane.normal.dot(b);
  CHECK(asgndist != bsgndist);
  return (a*bsgndist - b*asgndist)/(bsgndist - asgndist);
}

}              // anonymous methods

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polygon.
//------------------------------------------------------------------------------
std::string
polygon2string(const Polygon& poly) {
  ostringstream s;

  // Numbers of vertices.
  const auto nverts = poly.size();
  const auto nactive = count_if(poly.begin(), poly.end(),
                                [](const Vertex2d& x) { return x.comp >= 0; });
  set<int> usedVertices;

  // First dump the raw vertex info.
  s << "{\n";
  for (auto i = 0; i < nverts; ++i) {
    s << "  " << i << " " << poly[i].position
      << " [" << poly[i].neighbors.first << " " << poly[i].neighbors.second << "]\n";
  }
  s << "}\n";

  // // Go until we hit all the active vertices.
  // s << "[";
  // while (usedVertices.size() < nactive) {
  //   s << "[";

  //   // Look for the first active unused vertex.
  //   auto vstart = 0;
  //   while (vstart < nverts and
  //          (poly[vstart].comp < 0 or usedVertices.find(vstart) != usedVertices.end())) vstart++;
  //   CHECK(vstart < nverts);
  //   auto vnext = vstart;

  //   // Read out this loop.
  //   auto force = true;
  //   while (force or vnext != vstart) {
  //     s << " " << poly[vnext].position;
  //     force = false;
  //     usedVertices.insert(vnext);
  //     vnext = poly[vnext].neighbors.second;
  //   }
  //   s << "]";
  // }
  // s << "]";
  return s.str();
}

//------------------------------------------------------------------------------
// Convert Spheral::GeomPolygon -> PolyClipper::Polygon.
//------------------------------------------------------------------------------
void convertToPolygon(Polygon& polygon,
                      const Spheral::Dim<2>::FacetedVolume& Spheral_polygon) {
  TIME_PC2d_convertto.start();

  // Construct the vertices without connectivity first.
  const auto& coords = Spheral_polygon.vertices();
  const auto  nverts0 = coords.size();
  polygon.resize(nverts0);
  for (auto i = 0; i < nverts0; ++i) polygon[i].position = coords[i];

  // Build the connectivity.
  const auto& facets = Spheral_polygon.facets();
  int v1, v2, ivert1, ivert2;
  for (const auto& facet: facets) {
    v1 = facet.ipoint1();
    v2 = facet.ipoint2();
    CHECK(v1 < nverts0);
    CHECK(v2 < nverts0);
    polygon[v1].neighbors.second = v2;
    polygon[v2].neighbors.first  = v1;
  }

  CHECK(polygon.size() == nverts0);
  TIME_PC2d_convertto.stop();
}

//------------------------------------------------------------------------------
// Convert PolyClipper::Polygon -> Spheral::GeomPolygon.
//------------------------------------------------------------------------------
void convertFromPolygon(Spheral::Dim<2>::FacetedVolume& Spheral_polygon,
                        const Polygon& polygon) {
  TIME_PC2d_convertfrom.start();

  // Useful types.
  typedef Spheral::Dim<2>::FacetedVolume FacetedVolume;
  typedef Spheral::Dim<2>::Vector Vector;

  if (polygon.empty()) {

    Spheral_polygon = FacetedVolume();

  } else {

    // Numbers of vertices.
    const auto nverts = polygon.size();
    const auto nactive = count_if(polygon.begin(), polygon.end(),
                                  [](const Vertex2d& x) { return x.comp >= 0; });
    set<int> usedVertices;

    // Go until we hit all the active vertices.
    vector<Vector> coords(nactive);
    vector<vector<unsigned>> facets(nactive, vector<unsigned>(2));
    auto k = 0, loopStart = 0;
    while (usedVertices.size() < nactive) {

      // Look for the first active unused vertex.
      auto vstart = 0;
      while (vstart < nverts and
             (polygon[vstart].comp < 0 or usedVertices.find(vstart) != usedVertices.end())) vstart++;
      CHECK(vstart < nverts);
      auto vnext = vstart;

      // Read out this loop.
      auto force = true;
      while (force or vnext != vstart) {
        CHECK2(k < nactive, polygon2string(polygon));
        coords[k] = polygon[vnext].position;
        facets[k][0] = k;
        facets[k][1] = k + 1;
        ++k;
        force = false;
        usedVertices.insert(vnext);
        vnext = polygon[vnext].neighbors.second;
      }
      facets[k-1][1] = loopStart;
      loopStart = k;
    }
    CHECK(k == nactive);

    Spheral_polygon = FacetedVolume(coords, facets);

  }
  TIME_PC2d_convertfrom.stop();
}

//------------------------------------------------------------------------------
// Compute the zeroth and first moment of a Polygon.
//------------------------------------------------------------------------------
void moments(double& zerothMoment, Spheral::Dim<2>::Vector& firstMoment,
             const Polygon& polygon) {
  TIME_PC2d_moments.start();

  // Useful types.
  typedef Spheral::Dim<2>::Vector Vector;

  // Clear the result for accumulation.
  zerothMoment = 0.0;
  firstMoment = Vector::zero;

  // Walk the polygon, and add up our results triangle by triangle.
  if (not polygon.empty()) {
    const auto nverts = polygon.size();
    auto vfirst = 0;
    auto vprev = vfirst;
    auto vnext = vfirst;
    for (auto k = 0; k < nverts; ++k) {
      vnext = polygon[vnext].neighbors.second;
      const auto triA = (polygon[vprev].position - polygon[vfirst].position).cross(polygon[vnext].position - polygon[vfirst].position).z();
      zerothMoment += triA;
      firstMoment += triA * (polygon[vprev].position + polygon[vnext].position);
      vprev = vnext;
    }
    CHECK(zerothMoment != 0.0);
    firstMoment = firstMoment/(3.0*zerothMoment) + polygon[vfirst].position;
    zerothMoment *= 0.5;
  }
  TIME_PC2d_moments.stop();
}

//------------------------------------------------------------------------------
// Clip a polygon by planes.
//------------------------------------------------------------------------------
void clipPolygon(Polygon& polygon,
                 const std::vector<Plane2d>& planes) {
  TIME_PC2d_clip.start();

  // Useful types.
  typedef Spheral::Dim<2>::Vector Vector;

  // cerr << "Initial polygon: " << polygon2string(polygon) << endl;

  // Find the bounding box of the polygon.
  auto xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
  auto ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
  for (auto& v: polygon) {
    xmin = std::min(xmin, v.position[0]);
    xmax = std::max(xmax, v.position[0]);
    ymin = std::min(ymin, v.position[1]);
    ymax = std::max(ymax, v.position[1]);
  }

  // Loop over the planes.
  TIME_PC2d_planes.start();
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not polygon.empty()) {
    const auto& plane = planes[kplane++];
    // cerr << "Clip plane: " << plane.dist << " " << plane.normal << endl;

    // First check against the bounding box.
    auto boxcomp = compare(plane, xmin, ymin, xmax, ymax);
    auto above = boxcomp ==  1;
    auto below = boxcomp == -1;
    CHECK(not (above and below));

    // Check the current set of vertices against this plane.
    TIME_PC2d_checkverts.start();
    if (not (above or below)) {
      for (auto& v: polygon) {
        v.comp = compare(plane, v.position);
        if (v.comp == 1) {
          below = false;
        } else if (v.comp == -1) {
          above = false;
        }
      }
      CHECK(not (above and below));
    }
    TIME_PC2d_checkverts.stop();

    // Did we get a simple case?
    if (below) {
      // The polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes -- we're done.
      polygon.clear();

    } else if (not above) {

      // This plane passes through the polygon.
      // Insert any new vertices.
      TIME_PC2d_insertverts.start();
      vector<int> hangingVertices;
      int vprev, vnext, vnew;
      const auto nverts0 = polygon.size();
      for (auto v = 0; v < nverts0; ++v) {
        std::tie(vprev, vnext) = polygon[v].neighbors;

        if ((polygon[v].comp)*(polygon[vnext].comp) == -1) {
          // This pair straddles the plane and creates a new vertex.
          vnew = polygon.size();
          polygon.push_back(Vertex2d(segmentPlaneIntersection(polygon[v].position,
                                                              polygon[vnext].position,
                                                              plane),
                                     2));         // 2 indicates new vertex
          polygon[vnew].neighbors = {v, vnext};
          polygon[v].neighbors.second = vnew;
          polygon[vnext].neighbors.first = vnew;
          hangingVertices.push_back(vnew);
          // cerr << " --> Inserting new vertex @ " << polygon.back().position << endl;

        } else if (polygon[v].comp == 0 and 
                   (polygon[vprev].comp == -1 xor polygon[vnext].comp == -1)) {
          // This vertex is exactly in-plane, but has exactly one neighbor edge that will be entirely clipped.
          // No new vertex, but vptr will be hanging.
          hangingVertices.push_back(v);
          // cerr << " --> Hanging vertex @ " << polygon[v].position << endl;

        }
      }
      TIME_PC2d_insertverts.stop();

      // cerr << "After insertion: " << polygon2string(polygon) << endl;

      // For each hanging vertex, link to the neighbors that survive the clipping.
      // If there are more than two hanging vertices, we've clipped a non-convex face and need to check
      // how to hook up each section, possibly resulting in new faces.
      TIME_PC2d_hanging.start();
      CHECK(hangingVertices.size() % 2 == 0);
      if (true) { //(hangingVertices.size() > 2) {

        // Yep, more than one new edge here.
        const Vector direction(plane.normal.y(), -(plane.normal.x()));
        sort(hangingVertices.begin(), hangingVertices.end(), 
             [&](const int a, const int b) { return (polygon[a].position).dot(direction) < (polygon[b].position).dot(direction); });

        // Now the ordered pairs of these new vertices form the new edges.
        int v1, v2;
        for (auto k = 0; k < hangingVertices.size(); k += 2) {
          v1 = hangingVertices[k];
          v2 = hangingVertices[k + 1];
          polygon[v1].neighbors.second = v2;
          polygon[v2].neighbors.first  = v1;
        }

      } else {

        // Just hook across the vertices and we're done.
        for (auto v: hangingVertices) {
          std::tie(vprev, vnext) = polygon[v].neighbors;
          CHECK(polygon[v].comp == 0 or polygon[v].comp == 2);
          CHECK(polygon[vprev].comp == -1 xor polygon[vnext].comp == -1);

          if (polygon[vprev].comp == -1) {
            // We have to search backwards.
            while (polygon[vprev].comp == -1) {
              vprev = polygon[vprev].neighbors.first;
            }
            CHECK(vprev != v);
            polygon[v].neighbors.first = vprev;

          } else {
            // We have to search forward.
            while (polygon[vnext].comp == -1) {
              vnext = polygon[vnext].neighbors.second;
            }
            CHECK(vnext != v);
            polygon[v].neighbors.second = vnext;

          }
        }

      }

      // Remove the clipped vertices, compressing the polygon.
      TIME_PC2d_compress.start();

      // First, number the active vertices sequentially.
      xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
      ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
      auto i = 0;
      auto nkill = 0;
      for (auto& v: polygon) {
        if (v.comp < 0) {
          nkill++;
        } else {
          v.ID = i++;
          xmin = std::min(xmin, v.position[0]);
          xmax = std::max(xmax, v.position[0]);
          ymin = std::min(ymin, v.position[1]);
          ymax = std::max(ymax, v.position[1]);
        }
      }

      // Find the vertices to remove, and renumber the neighbors.
      if (nkill > 0) {
        vector<int> verts2kill;
        for (auto k = 0; k < polygon.size(); ++k) {
          if (polygon[k].comp < 0) {
            verts2kill.push_back(k);
          } else {
            polygon[k].neighbors.first = polygon[polygon[k].neighbors.first].ID;
            polygon[k].neighbors.second = polygon[polygon[k].neighbors.second].ID;
          }
        }
        Spheral::removeElements(polygon, verts2kill);
      }
      TIME_PC2d_compress.stop();

      // cerr << "After compression: " << polygon2string(polygon) << endl;

      // Is the polygon gone?
      if (polygon.size() < 3) polygon.clear();
      TIME_PC2d_hanging.stop();
    }
  }
  TIME_PC2d_planes.stop();
}

}
