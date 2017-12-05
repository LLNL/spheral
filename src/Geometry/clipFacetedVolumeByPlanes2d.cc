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

#include <iostream>
#include <iterator>

namespace Spheral {

using namespace std;

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
    if (vertexMask[i] == 1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
  }
  s << "\n"
    << "Clipped vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == -1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
  }
  s << "\n"
    << "Inactive vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == 0) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
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
void clipFacetedVolumeByPlanes(GeomPolygon& poly, 
                               const std::vector<GeomPlane<Dim<2>>>& planes) {

  // Useful types.
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::Vector Vector;
  typedef std::pair<int, int> Edge;       // Edges are pairs of vertex indices.  Edges are always built with (minVertID, maxVertID).
  typedef std::vector<int> Face;          // Faces are loops of edges.  We use the 1's complement to indicate an edge should be reversed 

  // Convert the polygon to a set of edge loops.
  auto vertices = poly.vertices();             // Note this is a copy!
  vector<Edge> edges;                          // Edges as pairs of vertex indices
  Face         face;                           // The rings of edges making up the polygon.
  vector<int> vertexMask(vertices.size(), 1);  // Mask to flag active/inactive vertices: 0->inactive, 1->active, -1->clip
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
      vector<int> newEdges;               // Any new edges we create.
      Face newface;                       // The newly clipped face as a set of edges.

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

          // Now edit our old edge in-place.
          edge = make_edge(vertID, v1);
          newface.push_back(edge.first == vertID ? iedge : ~iedge);

        } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {
          // v1 is clipped
          int vertID;
          insertVertex(vertices, vertID, vertexMask, v0, v1, p0, phat);

          // Now edit the old edge in place.
          edge = make_edge(v0, vertID);
          newface.push_back(edge.first == v0 ? iedge : ~iedge);

        } else {
          // Both v0 and v1 are clipped
          CHECK(vertexMask[v0] == -1 and vertexMask[v1] == -1);
        }
      }
      CHECK(vertices.size() == vertexMask.size());

      // We have walked and clipped the existing edges of the face.  First check if the face survived at all.
      if (not newface.empty()) {

        // OK, all the original edges of the face have been clipped, but we need to check for gaps along the clipping plane
        // that need to be connected with new edges.
        // cerr << "BLAGO!  ";
        // for (auto kedge = 0; kedge < newface.size(); ++kedge) cerr << " (" << startVertex(newface[kedge], edges) << " " << endVertex(newface[kedge], edges) << ")";
        // cerr << endl;
        for (auto kedge = 0; kedge < newface.size(); ++kedge) {
          const auto kedge1 = (kedge + 1) % newface.size();
          const auto v0 = endVertex(newface[kedge], edges);
          const auto v1 = startVertex(newface[kedge1], edges);
          if (v0 != v1) {
            // Yep, there's a gap here so insert a new edge.
            const int edgeID = edges.size();
            edges.push_back(make_edge(v0, v1));
            newEdges.push_back(                      (edges.back().first == v0 ? ~edgeID :  edgeID));
            newface.insert(newface.begin() + kedge1, (edges.back().first == v0 ?  edgeID : ~edgeID));
          }
        }
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
