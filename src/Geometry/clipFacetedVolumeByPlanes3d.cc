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
#include "Utilities/removeElements.hh"
#include "Utilities/DBC.hh"

#include <iostream>
#include <iterator>

namespace Spheral {

using namespace std;

namespace {

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
Dim<3>::Vector
segmentPlaneIntersection(const Dim<3>::Vector& a,       // line-segment begin
                         const Dim<3>::Vector& b,       // line-segment end
                         const Dim<3>::Vector& p,       // point in plane
                         const Dim<3>::Vector& phat) {  // plane unit normal

  const auto ab = b - a;
  const auto abhat = ab.unitVector();
  CHECK2(std::abs(abhat.dot(phat)) > 0.0, (abhat.dot(phat)) << " " << a << " " << b << " " << abhat << " " << phat);
  const double s = std::max(0.0, std::min(ab.magnitude(), (p - a).dot(phat)/(abhat.dot(phat))));
  CHECK2(s >= 0.0 and s <= ab.magnitude(), s << " " << ab.magnitude());
  const auto result = a + s*abhat;
  CHECK2(fuzzyEqual((result - p).dot(phat), 0.0, 1.0e-10),
         a << " " << b << " " << s << " " << result << " " << (result - p).dot(phat));
  return result;
}

//------------------------------------------------------------------------------
// Decide if the given edge plane intersection should give us a new vertex.
//------------------------------------------------------------------------------
inline
bool
insertVertex(std::vector<Dim<3>::Vector>& vertices,
             int& vertID,
             std::vector<int>& vertexMask,
             const int v0,
             const int v1,
             const Dim<3>::Vector& p0,
             const Dim<3>::Vector& phat) {

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
// Compare a plane and point (our built-in plane one has some issues).
//------------------------------------------------------------------------------
inline
int compare(const Dim<3>::Vector& planePoint,
            const Dim<3>::Vector& planeNormal,
            const Dim<3>::Vector& point) {
  const auto sgndist = planeNormal.dot(point - planePoint);
  if (std::abs(sgndist) < 1.0e-8) return 0;
  return sgn0(sgndist);
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
// When we edit a face in place we may flip the orientation required -- this 
// method is used to fix that orientation in the other face using the edge.
//------------------------------------------------------------------------------
inline
void fixOtherFaceEdgeOrientation(vector<vector<int>>& faces,
                                 const vector<pair<int, int>>& edgeFaces,
                                 const int faceID,
                                 const int edgeID,
                                 const int newEdgeOrientation) {
  REQUIRE(faceID >= 0 and faceID < faces.size());
  REQUIRE(edgeID >= 0 and edgeID < edgeFaces.size());
  REQUIRE(edgeFaces[edgeID].first == faceID or edgeFaces[edgeID].second == faceID);
  auto otherFaceID = (edgeFaces[edgeID].first == faceID ? edgeFaces[edgeID].second : edgeFaces[edgeID].first);
  // CHECK(otherFaceID > faceID);   // We should be the first to visit this edge.
  auto k = 0;
  while (k < faces[otherFaceID].size() and posID(faces[otherFaceID][k]) != edgeID) ++k;
  CHECK(k < faces[otherFaceID].size());
  faces[otherFaceID][k] = newEdgeOrientation;
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
// Return a formatted string with the current polyhedron.
//------------------------------------------------------------------------------
inline
string
poly2string(const vector<Dim<3>::Vector>& vertices,
            const vector<int>& vertexMask,
            const vector<pair<int, int>>& edges,
            const vector<vector<int>>& faces) {
  typedef Dim<3>::Vector Vector;
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
  for (auto iface = 0; iface < faces.size(); ++iface) {
    s << "Face " << iface << " : ";
    for (auto kedge = 0; kedge < faces[iface].size(); ++kedge) {
      const auto iedge = posID(faces[iface][kedge]);
      s <<  " ([" << faces[iface][kedge] << "->" << iedge << "] ";
      if (faces[iface][kedge] < 0) {
        s << edges[iedge].second << " " << edges[iedge].first << ")";;
      } else {
        s << edges[iedge].first << " " << edges[iedge].second << ")";
      }
    }
    s << "\n";
  }
  return s.str();
}

}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
void clipFacetedVolumeByPlanes(GeomPolyhedron& poly, 
                               const std::vector<GeomPlane<Dim<3>>>& planes) {

  // Useful types.
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::Vector Vector;
  typedef std::pair<int, int> Edge;       // Edges are pairs of vertex indices.  Edges are always built with (minVertID, maxVertID).
  typedef std::vector<int> Face;          // Faces are loops of edges.  We use the 1's complement to indicate an edge should be reversed in this face.

  // Convert the polyhedron to faces consisting of edge loops.
  auto vertices = poly.vertices();             // Note this is a copy!
  vector<Edge> edges;                          // Edges as pairs of vertex indices
  vector<Face> faces;                          // Faces that make up the polyhedron (loops of edges).
  vector<pair<int, int>> edgeFaces;            // The faces each edge is in (there should only be 2, hence the pair).
  vector<int> vertexMask(vertices.size(), 1);  // Mask to flag active/inactive vertices: 0->inactive, 1->active, -1->clip
  {
    const auto& facets = poly.facets();
    int iedge, iface = 0;
    for (const auto& facet: facets) {
      const auto& ipoints = facet.ipoints();
      const auto  npoints = ipoints.size();
      Face face;
      for (auto k = 0; k < npoints; ++k) {
        const auto edge = make_edge(ipoints[k], ipoints[(k + 1) % npoints]);
        auto itr = find(edges.begin(), edges.end(), edge);
        if (itr == edges.end()) {
          edges.push_back(edge);
          iedge = edges.size() - 1;
          edgeFaces.push_back(make_pair(iface, -1));
        } else {
          iedge = std::distance(edges.begin(), itr);
          CHECK(edgeFaces[iedge].first >= 0 and edgeFaces[iedge].second == -1);
          edgeFaces[iedge].second = iface;
        }
        face.push_back(edge.first == ipoints[k] ? iedge : ~iedge);
      }
      faces.push_back(face);
      ++iface;
    }
  }

  // // BLAGO
  // cerr << "----------------------------------------------------------------------" << endl
  //      << "Starting polyhedron: " << endl << poly2string(vertices, vertexMask, edges, faces) << endl;
  // // BLAGO

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not faces.empty()) {
    const auto& plane = planes[kplane++];
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();
    // cerr << "................................................................................" << endl
    //      << "Plane " << p0 << " " << phat << endl;

    // Check if the polyhedron is entirely clear of the plane (above or below).
    auto above = true;
    auto below = true;
    {
      auto k = 0;
      while (k < vertices.size()) {
        if (vertexMask[k] != 0) {
          const auto vcomp = compare(p0, phat, vertices[k]);
          if (vcomp >= 0) {
            below = false;
            vertexMask[k] = 1;          // Mark this vertex as keep.
          } else if (vcomp == -1) {
            above = false;
            vertexMask[k] = -1;         // Mark this vertex as clip.
          }
        }
        ++k;
      }
    }
    CHECK(not (above and below));

    // Did we get a simple case?
    if (below) {
      // The current polyhedron is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      edges.clear();
      faces.clear();

    } else if (not above) {
      vector<int> newEdges;                         // Any new edges we create.

      // This plane passes somewhere through the polyhedron, so we need to clip it.
      // Walk and edit each current face -- we'll handle adding the new faces after this pass.
      for (auto kface = 0; kface < faces.size(); ++kface) {
        // cerr << "Face " << kface << endl;
        auto& face = faces[kface];
        Face newface;

        // Walk the loop of edges in this face.
        for (auto kedge = 0; kedge < face.size(); ++kedge) {
          auto& edge = edges[posID(face[kedge])];
          auto  v0 = edge.first;
          auto  v1 = edge.second;
          if (face[kedge] < 0) std::swap(v0, v1);
          CHECK(vertexMask[v0] != 0 and vertexMask[v1] != 0);

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
            const int iedge = posID(face[kedge]);
            newface.push_back(edge.first == vertID ? iedge : ~iedge);

            // Since we edited an edge we have to double-check it's orientation in the other face using it.
            fixOtherFaceEdgeOrientation(faces, edgeFaces, kface, posID(newface.back()), ~newface.back());

          } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {
            // v1 is clipped
            int vertID;
            insertVertex(vertices, vertID, vertexMask, v0, v1, p0, phat);

            // Now edit the old edge in place.
            edge = make_edge(v0, vertID);
            const int iedge = posID(face[kedge]);
            newface.push_back(edge.first == v0 ? iedge : ~iedge);

            // Since we edited an edge we have to double-check it's orientation in the other face using it.
            fixOtherFaceEdgeOrientation(faces, edgeFaces, kface, posID(newface.back()), ~newface.back());

          } else {
            // Both v0 and v1 are clipped
            CHECK(vertexMask[v0] == -1 and vertexMask[v1] == -1);
          }
        }
        CHECK(vertices.size() == vertexMask.size());
        CHECK(edges.size() == edgeFaces.size());

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
              edgeFaces.push_back(make_pair(kface, -1));
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

        // The new face is complete, so replace the old one.
        face = newface;
        // cerr << poly2string(vertices, vertexMask, edges, faces) << endl;
      }

      // // BLAGO
      // {
      //   cerr << "newEdges: ";
      //   std::copy(newEdges.begin(), newEdges.end(), std::ostream_iterator<int>(std::cerr, " "));
      //   cerr << endl;
      // }
      // // BLAGO

      // Well, now all the starting faces of the polyhedron have been clipped by the plane, but we need
      // to cap off new edges created by the plane as new faces of the polyhedron.
      if (not newEdges.empty()) {

        // We start by ordering the new edges such that they form topological loops.  Note
        // there may be more than one loop if we are cutting a non-convex polyhedron.  We also know the directions
        // we want to walk these loops based on our sign convention for the ID in newEdges.
        const auto nNewEdges = newEdges.size();
        CHECK(nNewEdges >= 3);
        int v1;
        for (auto iedge = 0; iedge < nNewEdges - 1; ++iedge) {
          v1 = endVertex(newEdges[iedge], edges);
          auto jedge = iedge + 1;
          while (jedge < nNewEdges and v1 != startVertex(newEdges[jedge], edges)) ++jedge;
          if (jedge < nNewEdges) std::swap(newEdges[iedge + 1], newEdges[jedge]);
        }
        // cerr << "Ordered edges for capping: ";
        // for (auto iedge = 0; iedge != nNewEdges; ++iedge) {
        //   if (newEdges[iedge] >= 0) {
        //     cerr << " (" << edges[newEdges[iedge]].first << " " << edges[newEdges[iedge]].second << ")";
        //   } else {
        //     cerr << " (" << edges[~newEdges[iedge]].second << " " << edges[~newEdges[iedge]].first << ")";
        //   }
        // }
        // cerr << endl;

        // Now we can read out each loop from the ordered newEdges to make our new faces.
        auto iedge = 0;
        while (iedge < nNewEdges) {
          auto newFaceID = faces.size();
          Face newface(1U, newEdges[iedge]);
          v1 = endVertex(newEdges[iedge], edges);
          auto jedge = iedge + 1;
          while (jedge < nNewEdges and startVertex(newEdges[jedge], edges) == v1) {
            newface.push_back(newEdges[jedge]);
            v1 = endVertex(newEdges[jedge], edges);
            ++jedge;
          }
          CHECK(newface.size() >= 3);
          // if (!(startVertex(newface.front(), edges) == endVertex(newface.back(), edges))) {
          //   cerr << "Ordered edges for capping: ";
          //   for (auto iedge = 0; iedge != nNewEdges; ++iedge) {
          //     if (newEdges[iedge] >= 0) {
          //       cerr << " (" << edges[newEdges[iedge]].first << " " << edges[newEdges[iedge]].second << ")";
          //     } else {
          //       cerr << " (" << edges[~newEdges[iedge]].second << " " << edges[~newEdges[iedge]].first << ")";
          //     }
          //   }
          //   cerr << endl;
          //   cerr << "Bad new face: ";
          //   for (const auto iedge: newface) {
          //     const auto v0 = startVertex(iedge, edges);
          //     const auto v1 = endVertex(iedge, edges);
          //     cerr << " ([" << v0 << "] " << vertices[v0] << " -> [" << v1 << "] " << vertices[v1] << ")";
          //   }
          //   cerr << endl;
          // }
          CHECK(startVertex(newface.front(), edges) == endVertex(newface.back(), edges));
          faces.push_back(newface);
          iedge = jedge;

          // List the new face with all the new edges in it.
          for (const auto e: newface) {
            const auto eid = posID(e);
            CHECK(edgeFaces[eid].second == -1);
            edgeFaces[eid].second = newFaceID;
          }
        }
      }
    }

    // Deactivate any clipped vertices.
    for (auto k = 0; k < vertexMask.size(); ++k) vertexMask[k] = std::max(0, vertexMask[k]);
  }

  // If we have an empty polyhedron just finish it here.
  if (faces.empty()) {
    poly = GeomPolyhedron();
    return;
  }

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
    }
  }

  // Build the facet info.
  vector<vector<unsigned>> facets;
  for (const auto& face: faces) {
    if (not face.empty()) {
      facets.push_back(vector<unsigned>());
      for (const auto edgeID: face) {
        facets.back().push_back(old2new[startVertex(edgeID, edges)]);
      }
      CHECK(facets.back().size() == face.size());
      CHECK(facets.back().size() >= 3);
      CHECK(*max_element(facets.back().begin(), facets.back().end()) < newvertices.size());
    }
  }
  CHECK(facets.size() >= 4);

  // Now we can rebuild the polyhedron.
  poly = GeomPolyhedron(newvertices, facets);
}

}
