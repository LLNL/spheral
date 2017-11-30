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

namespace Spheral {

using namespace std;

namespace {

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
segmentPlaneIntersection(const Dim<3>::Vector& a,       // line-segment begin
                         const Dim<3>::Vector& b,       // line-segment end
                         const Dim<3>::Vector& p,       // point in plane
                         const Dim<3>::Vector& phat) {  // plane unit normal

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
int compare(const Dim<3>::Vector& planePoint,
            const Dim<3>::Vector& planeNormal,
            const Dim<3>::Vector& point) {
  return sgn0(planeNormal.dot(point - planePoint));
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
  vector<int> vertexMask(vertices.size(), 1);  // Mask to flag active/inactive vertices: 0->inactive, 1->active, -1->clip
  {
    const auto& facets = poly.facets();
    int iedge;
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
        } else {
          iedge = std::distance(edges.begin(), itr);
        }
        face.push_back(edge.first == ipoints[k] ? iedge : ~iedge);
      }
      faces.push_back(face);
    }
  }

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not faces.empty()) {
    const auto& plane = planes[kplane++];
    const auto& p0 = plane.point();
    const auto& phat = plane.normal();

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
      vector<int> faces2kill;                       // Any faces to be entirely clipped away.
      vector<int> newEdges;                         // Any new edges we create.
      vector<int> newEdgeVertex(edges.size(), -1);  // If an edge was clipped, the index of the newly created vertex.  -1 if the edge has no new vertices.

      // This plane passes somewhere through the polyhedron, so we need to clip it.
      // Walk and edit each current face -- we'll handle adding the new faces after this pass.
      for (auto kface = 0; kface < faces.size(); ++kface) {
        auto& face = faces[kface];
        Face newface;

        // Walk the loop of edges in this face.
        auto lastVertex = -1;
        auto hangingVertex = -1;
        for (auto kedge = 0; kedge < face.size(); ++kedge) {
          auto& edge = edges[posID(face[kedge])];
          auto  v0 = edge.first;
          auto  v1 = edge.second;
          if (face[kedge] < 0) std::swap(v0, v1);
          CHECK(vertexMask[v0] != 0 and vertexMask[v1] != 0);

          // Check the edge against the clip plane based on its vertices.
          if (vertexMask[v0] == 1 and vertexMask[v1] == 1) {
            // The edge got through unscathed.  Add it back to the face.
            newface.push_back(face[kedge]);
            lastVertex = v1;

            // One wrinkle -- if this edge was edited by a prior face loop, we could have an initial hanging node.
            if (kedge == 0 and newEdgeVertex[posID(face[kedge])] == v0) hangingVertex = v0;

          } else if (vertexMask[v0] == -1 and vertexMask[v1] == 1) {
            // v0 is clipped
            const auto newVertexID = vertices.size();
            vertices.push_back(segmentPlaneIntersection(vertices[v0], vertices[v1], p0, phat));
            vertexMask.push_back(1);
            CHECK(vertices.size() == newVertexID and vertexMask.size() == newVertexID);
            vertexMask[v0] = 0;
            newEdgeVertex[posID(face[kedge])] = newVertexID;

            // Note, if v0 was clipped we're just re-entering the allowed volume.  
            // Before we edit this existing edge, we first build the new edge bridging the gap to the last vertex if possible.
            // However, if we have not yet acquired a valid lastVertex, we list this as a hanging vertex to be resolved
            // after we walk the entire ring.
            if (lastVertex >= 0) {
              int newEdgeID = edges.size();
              edges.push_back(make_edge(lastVertex, newVertexID));
              newEdgeVertex.push_back(-1);
              newEdges.push_back(edges.back().first == lastVertex ? ~newEdgeID : newEdgeID);
              newface.push_back (edges.back().first == lastVertex ? newEdgeID : ~newEdgeID);
            } else {
              CHECK(hangingVertex == -1);  // This should only happen at most once per face!
              hangingVertex = newVertexID;
            }

            // Now edit our old edge in-place.
            edge = make_edge(newVertexID, v1);
            newface.push_back(~face[kedge]);
            lastVertex = v1;

          } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {
            // v1 is clipped
            const auto newVertexID = vertices.size();
            vertices.push_back(segmentPlaneIntersection(vertices[v0], vertices[v1], p0, phat));
            vertexMask.push_back(1);
            CHECK(vertices.size() == newVertexID and vertexMask.size() == newVertexID);
            vertexMask[v1] = 0;
            newEdgeVertex[posID(face[kedge])] = newVertexID;
            edge = make_edge(v0, newVertexID);
            newface.push_back(face[kedge]);
            lastVertex = newVertexID;

          } else {
            // Both v0 and v1 are clipped
            CHECK(vertexMask[v0] == -1 and vertexMask[v1] == -1);
            vertexMask[v0] = 0;
            vertexMask[v1] = 0;
          }
        }
        CHECK(vertices.size() == vertexMask.size());

        // We have walked and clipped the existing edges of the face.  First check if the face survived at all.
        if (newface.empty()) {

          // This face was entirely clipped, so we need to mark it for deletion.
          faces2kill.push_back(kface);

        } else if (newface.size() < face.size()) {
          CHECK(newface.size() >= 2);

          // The face was clipped.  If we have an unresolved hanging node hook it to the last vertex
          // making a new edge to close the face ring.
          if (hangingVertex >= 0) {
            CHECK(lastVertex >= 0);
            int newEdgeID = edges.size();
            edges.push_back(make_edge(lastVertex, hangingVertex));
            newEdgeVertex.push_back(-1);
            newEdges.push_back(edges.back().first == lastVertex ? ~newEdgeID : newEdgeID);
            newface.push_back (edges.back().first == lastVertex ? newEdgeID : ~newEdgeID);
          }

          // The newly clipped face is complete, so replace the old one.
          face = newface;
        }
      }
      CHECK(newEdgeVertex.size() == edges.size());
      CHECK(std::count_if(newEdgeVertex.begin(), newEdgeVertex.end(), [](const int& x) { return x >= 0; }) == newEdges.size());

      // Gack any faces we have clipped.
      removeElements(faces, faces2kill);

      // Well, now all the starting faces of the polyhedron have been clipped by the plane, but we need
      // to cap off new edges created by the plane as new faces of the polyhedron.
      if (not newEdges.empty()) {

        // We start by ordering the new edges such that they form topological loops.  Note
        // there may be more than one loop if we are cutting a non-convex polyhedron.  We also know the directions
        // we want to walk these loops based on the sign of the ID in newEdges.
        const auto nNewEdges = newEdges.size();
        CHECK(nNewEdges >= 3);
        int v1;
        for (auto iedge = 0; iedge < nNewEdges - 1; ++iedge) {
          v1 = newEdges[iedge] >= 0 ? edges[newEdges[iedge]].second : edges[~newEdges[iedge]].first;
          auto jedge = iedge + 1;
          while (jedge < nNewEdges and 
                 (v1 != edges[posID(newEdges[jedge])].first and v1 != edges[posID(newEdges[jedge])].second)) ++jedge;
          if (jedge < nNewEdges) std::swap(newEdges[iedge + 1], newEdges[jedge]);
        }

        // Now we can read out each loop from the ordered newEdges to make our new faces.
        auto iedge = 0;
        while (iedge < nNewEdges) {
          Face newface(1U, newEdges[iedge]);
          v1 = newEdges[iedge] >= 0 ? edges[newEdges[iedge]].second : edges[~newEdges[iedge]].first;
          auto jedge = iedge + 1;
          while (jedge < nNewEdges and
                 (edges[posID(newEdges[jedge])].first == v1 or edges[posID(newEdges[jedge])].second == v1)) {
            newface.push_back(newEdges[jedge]);
            v1 = newEdges[jedge] >= 0 ? edges[newEdges[jedge]].second : edges[~newEdges[jedge]].first;
            ++jedge;
          }
          CHECK(newface.size() >= 3);
          faces.push_back(newface);
          iedge = jedge + 1;
        }
      }
    }
  }

  // If we have an empty polyhedron just finish it here.
  if (faces.empty()) {
    poly = GeomPolyhedron();
    return;
  }

  // Build a list of the vertices that are active, and a map of old->new vertex index.
  CHECK(*min_element(vertexMask.begin(), vertexMask.end()) == 0 and
        *max_element(vertexMask.begin(), vertexMask.end()) == 1);
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
  vector<vector<unsigned>> facets(faces.size());
  for (auto kface = 0; kface < faces.size(); ++kface) {
    const auto& face = faces[kface];
    for (const auto edgeID: face) {
      facets[kface].push_back(edgeID >= 0 ? old2new[edges[edgeID].first] : old2new[edges[~edgeID].second]);
    }
    CHECK(facets[kface].size() == face.size());
  }

  // Now we can rebuild the polyhedron.
  poly = GeomPolyhedron(newvertices, facets);
}

}
