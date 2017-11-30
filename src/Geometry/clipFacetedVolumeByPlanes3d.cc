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
  vector<int> vertexMask(vertices.size(), 1);  // Mask to flag active/inactive vertices.
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
      while (k < vertices.size() and (above or below)) {
        if (vertexMask[k] == 1) {
          const auto vcomp = compare(p0, phat, vertices[k]);
          if (vcomp >= 0) {
            below = false;
          } else if (vcomp == -1) {
            above = false;
            vertexMask[k] = -1;         // Mark this vertex as newly clipped.
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
      vector<unsigned> faces2kill;     // Any faces to be entirely clipped away.
      vector<unsigned> newedges;       // Any new edges we create.

      // This plane passes somewhere through the polyhedron, so we need to clip it.
      // Walk each current face.  Note we are not yet worried about adding the new faces.
      for (auto kface = 0; kface < faces.size(); ++kface) {
        auto& face = faces[kface];
        Face newface;

        // Walk the loop of edges in this face.
        for (auto kedge = 0; kedge < face.size(); ++kedge) {
          auto& edge = edges[posID(face[kedge])];
          auto  v0 = edge.first;
          auto  v1 = edge.second;
          CHECK(vertexMask[v0] != 0 and vertexMask[v0] != 0);

          // Is this edge clipped?
          if (vertexMask[v0] == 1 and vertexMask[v1] == 1) {
            // The edge got through unscathed.  Just add it back to the face.
            newface.push_back(face[kedge]);

          } else if (vertexMask[v0] == -1 and vertexMask[v1] == 1) {
            // v0 is clipped, but not v1.  Do the intersection and create a replacement edge.
            vertices.push_back(segmentPlaneIntersection(vertices[v0], vertices[v1], p0, phat));
            vertexMask.push_back(2);
            vertexMask[v0] = 0;
            edge = make_edge(v1, vertices.size() - 1);
            newface.push_back(~face[kedge]);

          } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {
            // v1 is clipped, but not v0.  Do the intersection and create a replacement edge.
            vertices.push_back(segmentPlaneIntersection(vertices[v0], vertices[v1], p0, phat));
            vertexMask.push_back(2);
            vertexMask[v1] = 0;
            edge = make_edge(v0, vertices.size() - 1);
            newface.push_back(face[kedge]);

          } else {
            // If none of the above branches hit, this edge is entirely clipped.
            CHECK(vertexMask[v0] == -1 and vertexMask[v1] == -1);
            vertexMask[v0] = 0;
            vertexMask[v1] = 0;
          }
        }
        CHECK(vertices.size() == vertexMask.size());

        // We have walked and clipped the existing edges of the face.  First check if the face survived at all.
        if (newface.size() == 0) {

          // This face was entirely clipped, so we need to mark it for deletion.
          faces2kill.push_back(kface);

        } else if (newface.size() < face.size()) {
          CHECK(newface.size() >= 2);

          // The face was clipped, so there must be some new hanging nodes we need to close with a new edge.
          // Since we're only clipping by a single plane in this pass there should only be a single gap.
          auto kedge = 0;
          while (kedge < newface.size() and 
                 not (vertexMask[edges[newface[posID(newface[kedge])]].first]  == 2 or
                      vertexMask[edges[newface[posID(newface[kedge])]].second] == 2)) ++kedge;
          CHECK(kedge < newface.size());
          
          // There is a gap in the new face between kedge and kedge+1, so create the new edge by connecting the hanging
          // new nodes.
          const auto& edge0 = edges[newface[posID(newface[kedge])]];
          const auto& edge1 = edges[newface[posID(newface[(kedge + 1) % newface.size()])]];
          CHECK(vertexMask[edge1.first] == 2 or vertexMask[edge1.second] == 2);
          const auto v0 = (vertexMask[edge0.first] == 2 ? edge0.first : edge0.second);
          const auto v1 = (vertexMask[edge1.first] == 2 ? edge1.first : edge1.second);
          CHECK(vertexMask[v0] == 2 and vertexMask[v1] == 2);    // Both new hanging nodes.
          edges.push_back(make_edge(v0, v1));
          newface.insert(newface.begin() + kedge + 1, edges.back().first == v0 ? (edges.size() - 1) : ~(edges.size() - 1));
          newedges.push_back(edges.size() - 1);

          // The hanging nodes have been hooked up, so flag them as active but not new.
          vertexMask[v0] = 1;
          vertexMask[v1] = 1;

          // OK, the newly clipped face is complete, so replace the old one.
          face = newface;
        }
      }

      // Gack any faces we have clipped.
      removeElements(faces, faces2kill);

      // Well, now all the starting faces of the polyhedron have been clipped by the plane, but we need
      // to cap off new edges created by the plane as new faces of the polyhedron.
      if (newedges.size() > 0) {
        vector<int> newEdgeMask(newedges.size(), 1);
        
        // Start with the first untouched new edge.
        auto startItr = std::find(newEdgeMask.begin(), newEdgeMask.end(), 1);
        if (startItr < newEdgeMask.end()) {

          // Build a new face, starting with this edge and walking connected edges until we close the loop.
          Face newface;

        }
      }
    }
  }

  // Build the final polyhedron, and we're done.
}

}
