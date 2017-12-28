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
#include "Utilities/Timer.hh"

#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>

// Declare the timers.
extern Timer TIME_CFV2d;
extern Timer TIME_CFV2d_convertfrom;
extern Timer TIME_CFV2d_planes;
extern Timer TIME_CFV2d_checkverts;
extern Timer TIME_CFV2d_insertverts;
extern Timer TIME_CFV2d_hanging;
extern Timer TIME_CFV2d_convertto;

namespace Spheral {

using namespace std;

namespace {

//------------------------------------------------------------------------------
// The vertex struct, which we use to encode the polygon.
//------------------------------------------------------------------------------
struct Vertex {
  typedef Dim<2>::Vector Vector;
  Vector position;
  std::pair<Vertex*, Vertex*> neighbors;
  int ID, comp;
  Vertex(): position(), neighbors(), ID(-1), comp(1) {}
  Vertex(const Vector& pos): position(pos), neighbors(), ID(-1), comp(1) {}
  Vertex(const Vector& pos, const int c): position(pos), neighbors(), ID(-1), comp(c) {}
};

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

}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
void clipFacetedVolumeByPlanes(GeomPolygon& poly0,
                               const std::vector<GeomPlane<Dim<2>>>& planes) {

  // Useful types.
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::Vector Vector;

  // Convert the input polygon to a loop of Vertex's.
  TIME_CFV2d.start();
  TIME_CFV2d_convertfrom.start();
  list<Vertex> poly;
  {
    const auto& vertPositions = poly0.vertices();
    const auto& facets = poly0.facets();
    const auto  nverts0 = vertPositions.size();
    vector<Vertex*> id2vert(nverts0, NULL);
    int v1, v2, ivert1, ivert2;
    for (const auto& facet: facets) {
      v1 = facet.ipoint1();
      v2 = facet.ipoint2();
      CHECK(v1 < nverts0);
      CHECK(v2 < nverts0);
      if (id2vert[v1] == NULL) {
        poly.push_back(Vertex(vertPositions[v1]));
        id2vert[v1] = &poly.back();
      }
      if (id2vert[v2] == NULL) {
        poly.push_back(Vertex(vertPositions[v2]));
        id2vert[v2] = &poly.back();
      }
      id2vert[v1]->neighbors.second = id2vert[v2];
      id2vert[v2]->neighbors.first =  id2vert[v1];
    }
    CHECK(poly.size() == nverts0);
  }

  // Initialize the set of active vertices.
  set<Vertex*> activeVertices;
  for (auto& v: poly) activeVertices.insert(&v);
  TIME_CFV2d_convertfrom.stop();

  // Loop over the planes.
  TIME_CFV2d_planes.start();
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not activeVertices.empty()) {
    const auto& plane = planes[kplane++];
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();
    // cerr << "Plane " << kplane << " : " << p0 << " " << phat << endl;

    // Check the current set of vertices against this plane.
    TIME_CFV2d_checkverts.start();
    auto above = true;
    auto below = true;
    for (auto vptr: activeVertices) {
      vptr->comp = compare(p0, phat, vptr->position);
      if (vptr->comp >= 0) {
        below = false;
      } else if (vptr->comp == -1) {
        above = false;
      }
    }
    CHECK(not (above and below));
    TIME_CFV2d_checkverts.stop();

    // Did we get a simple case?
    if (below) {
      // Polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      activeVertices.clear();

    } else if (not above) {

      // This plane passes through the polygon.
      // Insert any new vertices.
      TIME_CFV2d_insertverts.start();
      vector<Vertex*> hangingVertices, newVertices;
      Vertex *vprev, *vnext;
      for (auto vptr: activeVertices) {
        std::tie(vprev, vnext) = vptr->neighbors;

        if ((vptr->comp)*(vnext->comp) == -1) {
          // This pair straddles the plane and creates a new vertex.
          poly.push_back(Vertex(segmentPlaneIntersection(vptr->position,
                                                         vnext->position,
                                                         p0,
                                                         phat),
                                0));
          poly.back().neighbors.first = vptr;
          poly.back().neighbors.second = vnext;
          vptr->neighbors.second = &poly.back();
          vnext->neighbors.first = &poly.back();
          hangingVertices.push_back(&poly.back());
          newVertices.push_back(&poly.back());
          // cerr << "  --> Inserted vertex @ " << poly.back().position << endl;

        } else if (vptr->comp == 0 and 
                   (vprev->comp == -1 xor vnext->comp == -1)) {
          // This vertex is exactly in-plane, but has exactly one neighbor edge that will be entirely clipped.
          // No new vertex, but vptr will be hanging.
          hangingVertices.push_back(vptr);
          // cerr << "  --> Adding in-plane node @ " << vptr->position << endl;

        }
      }
      TIME_CFV2d_insertverts.stop();

      // For each hanging vertex, link to the neighbors that survive the clipping.
      TIME_CFV2d_hanging.start();
      for (auto vptr: hangingVertices) {
        std::tie(vprev, vnext) = vptr->neighbors;
        CHECK(vptr->comp == 0);
        CHECK(vprev->comp == -1 xor vnext->comp == -1);

        if (vprev->comp == -1) {
          // We have to search backwards.
          while (vprev->comp == -1) {
            activeVertices.erase(vprev);
            vprev = vprev->neighbors.first;
          }
          CHECK(vprev != vptr);
          vptr->neighbors.first = vprev;

        } else {
          // We have to search forward.
          while (vnext->comp == -1) {
            activeVertices.erase(vnext);
            vnext = vnext->neighbors.second;
          }
          CHECK(vnext != vptr);
          vptr->neighbors.second = vnext;

        }
      }

      // Add the new vertices to the active set.
      activeVertices.insert(newVertices.begin(), newVertices.end());
      CHECK(activeVertices.size() >= 3);
      TIME_CFV2d_hanging.stop();
    }
  }
  TIME_CFV2d_planes.stop();

  // Now rebuild the polygon, and we're done.
  TIME_CFV2d_convertto.start();
  if (activeVertices.empty()) {
    poly0 = GeomPolygon();
  } else {

    const auto nverts = activeVertices.size();
    CHECK(nverts >= 3);
    vector<Vector> coords(nverts);
    vector<vector<unsigned>> facets(nverts, vector<unsigned>(2));
    auto vptr = *activeVertices.begin();
    for (auto k = 0; k < nverts; ++k) {
      coords[k] = vptr->position;
      facets[k][0] = k;
      facets[k][1] = (k + 1) % nverts;
      vptr = vptr->neighbors.second;
    }
    CHECK(vptr == *activeVertices.begin());
    poly0 = GeomPolygon(coords, facets);
  }
  TIME_CFV2d_convertto.stop();

  TIME_CFV2d.stop();
}

}
