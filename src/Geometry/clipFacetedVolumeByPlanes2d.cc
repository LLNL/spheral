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

#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>

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

//------------------------------------------------------------------------------
// Build an edge as a unique std::pair.
//------------------------------------------------------------------------------
inline
std::pair<int, int>
make_edge(const int a, const int b) {
  return std::make_pair(std::min(a, b), std::max(a, b));
}

//------------------------------------------------------------------------------
// Return the ID for the given encoded value (1's complement crap).
//------------------------------------------------------------------------------
inline
int
posID(const int x) {
  return (x >= 0 ? x : ~x);
}

//------------------------------------------------------------------------------
// Grab the starting vertex for the edge corresponding to the encoded id.
//------------------------------------------------------------------------------
inline
int
startVertex(const int id, const vector<pair<int, int>>& edges) {
  if (id < 0) {
    CHECK(~id < edges.size());
    return edges[~id].second;
  } else {
    CHECK(id < edges.size());
    return edges[id].first;
  }
}

//------------------------------------------------------------------------------
// Grab the end vertex for the edge corresponding to the encoded id.
//------------------------------------------------------------------------------
inline
int
endVertex(const int id, const vector<pair<int, int>>& edges) {
  if (id < 0) {
    CHECK(~id < edges.size());
    return edges[~id].first;
  } else {
    CHECK(id < edges.size());
    return edges[id].second;
  }
}

//------------------------------------------------------------------------------
// Decide if the given edge plane intersection should give us a new vertex.
//------------------------------------------------------------------------------
inline
bool
insertVertex(std::vector<Dim<2>::Vector>& vertices,
             int& vertID,
             std::vector<int>& vertexMask,
             const int v0,
             const int v1,
             const Dim<2>::Vector& p0,
             const Dim<2>::Vector& phat) {

  // Where would we like the vertex?
  const auto v = segmentPlaneIntersection(vertices[v0], vertices[v1], p0, phat);

  // Is this degenerate with any existing active vertices?
  bool result;
  vertID = 0;
  while (vertID < vertices.size() and 
         not (vertexMask[vertID] == 1 and (vertices[vertID] - v).magnitude2() < 1.0e-16)) ++vertID;
  if (vertID == vertices.size()) {
    vertices.push_back(v);
    vertexMask.push_back(1);
    result = true;
  } else {
    result = false;
  }
  ENSURE(vertID < vertices.size());
  ENSURE(vertices.size() == vertexMask.size());
  ENSURE(vertexMask[vertID] == 1);
  return result;
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

//------------------------------------------------------------------------------
// Return a formatted string with the current polygon.
//------------------------------------------------------------------------------
inline
string
poly2string(const vector<Dim<2>::Vector>& vertices,
            const vector<int>& vertexMask,
            const vector<pair<int, int>>& edges,
            const vector<int>& face) {
  typedef Dim<2>::Vector Vector;
  ostringstream s;
  REQUIRE(vertices.size() == vertexMask.size());
  s << "Vertices: ";
  copy(vertices.begin(), vertices.end(), ostream_iterator<Vector>(s, " "));
  s << "\n"
    << "Vertex mask: ";
  copy(vertexMask.begin(), vertexMask.end(), ostream_iterator<int>(s, " "));
  s << "\n"
    << "Active vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == 1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << ")";
  }
  s << "\n"
    << "Clipped vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == -1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << ")";
  }
  s << "\n"
    << "Inactive vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == 0) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << ")";
  }
  s << "\n";
  s << "Face : ";
  for (auto kedge = 0; kedge < face.size(); ++kedge) {
    s <<  " ([" << face[kedge] << "->" << posID(face[kedge]) << "] ";
    s << startVertex(face[kedge], edges) << " " << endVertex(face[kedge], edges) << ")";
  }
  s << "\n";
  return s.str();
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

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not activeVertices.empty()) {
    const auto& plane = planes[kplane++];
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();

    // Check the current set of vertices against this plane.
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

    // Did we get a simple case?
    if (below) {
      // Polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      activeVertices.clear();

    } else if (not above) {

      // This plane passes through the polygon.
      // Insert any new vertices.
      vector<Vertex*> newVertices;
      for (auto vptr: activeVertices) {
        auto* vnext = vptr->neighbors.second;
        if ((vptr->comp)*(vnext->comp) == -1) {   // Does this pair straddle the plane?
          poly.push_back(Vertex(segmentPlaneIntersection(vptr->position,
                                                         vnext->position,
                                                         p0,
                                                         phat),
                                0));
          poly.back().neighbors.first = vptr;
          poly.back().neighbors.second = vnext;
          vptr->neighbors.second = &poly.back();
          vnext->neighbors.first = &poly.back();
          newVertices.push_back(&poly.back());
          activeVertices.insert(&poly.back());
        }
      }

      // For each new vertex, link to the neighbors that survive the clipping.
      for (auto vptr: newVertices) {
        CHECK(vptr->comp == 0);
        CHECK((vptr->neighbors.first->comp)*(vptr->neighbors.second->comp) == -1);
        if (vptr->neighbors.first->comp == -1) {
          // We have to search backwards.
          auto vprior = vptr->neighbors.first;
          while (vprior->comp == -1) {
            activeVertices.erase(vprior);
            vprior = vprior->neighbors.first;
          }
          CHECK(vprior != vptr);
          vptr->neighbors.first = vprior;
        } else {
          // We have to search forward.
          auto vnext = vptr->neighbors.second;
          while (vnext->comp == -1) {
            activeVertices.erase(vnext);
            vnext = vnext->neighbors.second;
          }
          CHECK(vnext != vptr);
          vptr->neighbors.second = vnext;
        }
      }
      CHECK(activeVertices.size() >= 3);
    }
  }

  // Now rebuild the polygon, and we're done.
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
}

}
