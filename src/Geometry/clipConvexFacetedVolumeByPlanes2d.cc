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
// This version is specialized for convex polygons/polyhedra.
//
// Created by J. Michael Owen, Wed Dec 13 15:28:09 PST 2017
//----------------------------------------------------------------------------//

#include "clipConvexFacetedVolumeByPlanes.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>
#include <iterator>
#include <algorithm>

namespace Spheral {


namespace {

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
void clipConvexFacetedVolumeByPlanes(GeomPolygon& poly, 
                                     const std::vector<GeomPlane<Dim<2>>>& planes) {

  // Useful types.
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::Vector Vector;
  typedef std::pair<int, int> Edge;       // Edges are pairs of vertex indices.  Edges are always built with (minVertID, maxVertID).
  typedef std::vector<int> Face;          // Faces are loops of edges.  We use the 1's complement to indicate an edge should be reversed 

  // We'd better be convex!
  VERIFY2(poly.convex(),
          "clipConvexFacetedVolumeByPlanes2d ERROR: require a convex polygon!");

  // Convert the polygon to a set of edge loops.
  auto vertices = poly.vertices();             // Note this is a copy!
  vector<Edge> edges;                          // edges as pairs of vertex indices
  Face         face;                           // the rings of edges making up the polygon.
  vector<int> vertexMask(vertices.size(), 1);  // mask to flag active/inactive vertices: 0->inactive, 1->active, -1->clip
  {
    const auto& facets = poly.facets();
    int iedge;
    for (const auto& facet: facets) {
      const auto edge = make_edge(facet.ipoint1(), facet.ipoint2());
      CHECK(find(edges.begin(), edges.end(), edge) == edges.end());
      iedge = edges.size();
      edges.push_back(edge);
      face.push_back(edge.first == facet.ipoint1() ? iedge : ~iedge);
    }
  }
  // cerr << "Initial polygon: "<< endl
  //      << poly2string(vertices, vertexMask, edges, face);

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not face.empty()) {
    const auto& plane = planes[kplane++];
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();

    // Check if the polygon is entirely clear of the plane (above or below).
    auto above = true;
    auto below = true;
    {
      auto k = 0;
      while (k < vertices.size()) {
        if (vertexMask[k] != 0) {
          const auto vcomp = compare(p0, phat, vertices[k]);
          if (vcomp >= 0) {
            below = false;
            vertexMask[k] = 1;           // Mark this vertex as keep.
          } else if (vcomp == -1) {
            above = false;
            vertexMask[k] = -1;          // Mark this vertex as clip.
          }
        }
        ++k;
      }
    }
    CHECK(not (above and below));

    // Did we get a simple case?
    if (below) {
      // Polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      edges.clear();
      face.clear();

    } else if (not above) {
      Face newface;                       // The newly clipped face as a set of edges.
      vector<int> newEdges;               // Any new edges we create.
      vector<pair<int, int>> newVertices; // And new vertices we create: store (vertex id, index of edge in the new face that was clipped).

      // The plane passes somewhere through the polygon, so we need to walk the ring and modify it.
      for (auto kedge = 0; kedge < face.size(); ++kedge) {
        const int iedge = posID(face[kedge]);
        auto& edge = edges[iedge];
        auto  v0 = startVertex(face[kedge], edges);
        auto  v1 = endVertex(face[kedge], edges);

        // Check the edge against the clip plane based on its vertices.
        if (vertexMask[v0] == 1 and vertexMask[v1] == 1) {
          // This edge got through unscathed -- add it back to the face.
          newface.push_back(face[kedge]);

        } else if (vertexMask[v0] == -1 and vertexMask[v1] == 1) {
          // v0 is clipped
          // Check if we're inserting a new vertex.
          int vertID;
          insertVertex(vertices, vertID, vertexMask, v0, v1, p0, phat);
          newVertices.push_back(make_pair(vertID, newface.size()));

          // Now edit our old edge in-place.
          edge = make_edge(vertID, v1);
          newface.push_back(edge.first == vertID ? iedge : ~iedge);

        } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {
          // v1 is clipped
          int vertID;
          insertVertex(vertices, vertID, vertexMask, v0, v1, p0, phat);
          newVertices.push_back(make_pair(vertID, newface.size()));

          // Now edit the old edge in place.
          edge = make_edge(v0, vertID);
          newface.push_back(edge.first == v0 ? iedge : ~iedge);

        } else {
          // Both v0 and v1 are clipped
          CHECK(vertexMask[v0] == -1 and vertexMask[v1] == -1);
        }
      }
      CHECK(vertices.size() == vertexMask.size());

      // Check if we clipped any edges.
      if (not newVertices.empty()) {
        CHECK(newVertices.size() == 2);

        // cerr << "Before adding new edges: "<< endl
        //      << poly2string(vertices, vertexMask, edges, newface);

        // OK, all the original edges of the face have been clipped, but there are gaps along the clipping
        // plane which require new edges be constructed.
        // Because this is a convex face, we know there can be at most two points.  Orient 'em correctly.
        const Vector direction(phat.y(), -phat.x());
        if ((vertices[newVertices[0].first] - p0).dot(direction) > (vertices[newVertices[1].first] - p0).dot(direction)) std::swap(newVertices[0], newVertices[1]);

        // Now the ordered pairs of these new vertices form the new edges.
        const int v0 = newVertices[0].first;
        const int iedge = edges.size();
        const auto newedge = make_edge(newVertices[0].first, newVertices[1].first);
        edges.push_back(newedge);
        newface.push_back(newedge.first == v0 ? iedge : ~iedge);
        CHECK(newface.size() >= 3);

        // Verify the new face forms a proper topological ring.
        BEGIN_CONTRACT_SCOPE
        {
          for (auto kedge = 0; kedge < newface.size(); ++kedge) {
            const auto kedge1 = (kedge + 1) % newface.size();
            CHECK(endVertex(newface[kedge], edges) == startVertex(newface[kedge1], edges));
          }
        }
        END_CONTRACT_SCOPE
      }

      // Deactivate any clipped vertices.
      for (auto k = 0; k < vertexMask.size(); ++k) vertexMask[k] = std::max(0, vertexMask[k]);

      // The new face is complete, so replace the old one.
      face = newface;
    }
  }
  // cerr << "Final clipped polygon: "<< endl
  //      << poly2string(vertices, vertexMask, edges, face);

  // Now rebuild the polygon, and we're done.
  if (face.empty()) {
    poly = GeomPolygon();
  } else {
    CHECK(face.size() >= 3);
    // cerr << "Building new polygon:" << endl
    //      << "    vertices: ";

    // Build a list of the vertices that are active, and a map of old->new vertex index.
    CHECK(vertexMask.empty() or 
          (*min_element(vertexMask.begin(), vertexMask.end()) >= 0 and
           *max_element(vertexMask.begin(), vertexMask.end()) <= 1));
    vector<Vector> newvertices;
    newvertices.reserve(vertices.size());
    vector<int> old2new(vertices.size());
    for (auto i = 0; i < vertices.size(); ++i) {
      if (vertexMask[i] == 1) {
        old2new[i] = newvertices.size();
        newvertices.push_back(vertices[i]);
        // cerr << " ([" << old2new[i] << "] " << newvertices.back() << ")";
      }
    }
    // cerr << endl
    //      << "    facets: ";
    // Build the facet info.
    vector<vector<unsigned>> facets(face.size(), vector<unsigned>(2));;
    for (auto kedge = 0; kedge < face.size(); ++kedge) {
      facets[kedge][0] = old2new[startVertex(face[kedge], edges)];
      facets[kedge][1] = old2new[endVertex(face[kedge], edges)];
      // cerr << "(" << facets[kedge][0] << " " << facets[kedge][1] << ")";
    }
    // for (const auto edgeID: face) {
    //   facets.push_back(vector<unsigned>(2));
    //   facets.back()[0] = old2new[startVertex(edgeID, edges)];
    //   facets.back()[1] = old2new[endVertex(edgeID, edges)];
    //   cerr << "(" << facets.back()[0] << " " << facets.back()[1] << ")";
    // }
    CHECK(facets.size() >= 3);
    // cerr << endl;

    // Now we can rebuild the polyhedron.
    poly = GeomPolygon(newvertices, facets);
  }
}

}
