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
#include "Utilities/timingUtilities.hh"

#include "easy/profiler.h"

#include <iostream>
#include <iterator>

namespace Spheral {

using namespace std;

namespace {

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
  return id < 0 ? edges[~id].second : edges[id].first;
  // if (id < 0) {
  //   CHECK(~id < edges.size());
  //   return edges[~id].second;
  // } else {
  //   CHECK(id < edges.size());
  //   return edges[id].first;
  // }
}

//------------------------------------------------------------------------------
// Grab the end vertex for the edge corresponding to the encoded id.
//------------------------------------------------------------------------------
inline
int
endVertex(const int id, const vector<pair<int, int>>& edges) {
  return id < 0 ? edges[~id].first : edges[id].second;
  // if (id < 0) {
  //   CHECK(~id < edges.size());
  //   return edges[~id].first;
  // } else {
  //   CHECK(id < edges.size());
  //   return edges[id].second;
  // }
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
// Sort a set of edges into topological loops.
//------------------------------------------------------------------------------
inline
void
sortInTopologicalRings(vector<int>& loops,
                       const vector<pair<int, int>>& edges) {
  const auto nedges = loops.size();
  CHECK(nedges >= 3);
  int v1, iedge, jedge;
  for (iedge = 0; iedge < nedges - 1; ++iedge) {
    v1 = endVertex(loops[iedge], edges);
    jedge = iedge + 1;
    while (jedge < nedges and v1 != startVertex(loops[jedge], edges)) ++jedge;
    if (jedge < nedges) std::swap(loops[iedge + 1], loops[jedge]);
  }

  // In some degenerate cases it's possible to have an edge we walk in reverse,
  // essentially negating it.  Check for this and remove any instances.
  vector<int> edges2kill;
  for (iedge = 0; iedge < nedges; ++iedge) {
    jedge = (iedge + 1) % nedges;
    if (loops[iedge] == ~loops[jedge]) {
      edges2kill.push_back(iedge);
      edges2kill.push_back(jedge);
    }
  }
  removeElements(loops, edges2kill);
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
    << "Above vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == 1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
  }
  s << "\n"
    << "Below vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == -1) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
  }
  s << "\n"
    << "In-plane vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == 0) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
  }
  s << "\n"
    << "Inactive vertices: ";
  for (auto i = 0; i < vertices.size(); ++i) {
    if (vertexMask[i] == -2) s << " ([" << i << "] " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << ")";
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

  // Timer.
  EASY_FUNCTION(profiler::colors::Magenta);

  // // The timing variables.
  // double tconvertfrom = 0.0,
  //         tvertexmask = 0.0,
  //           tedgemask = 0.0,
  //           tfaceclip = 0.0, 
  //   tfaceclip_edges = 0.0,
  //   tfaceclip_clear = 0.0,
  //   tfaceclip_dups = 0.0,
  //   tfaceclip_sortline = 0.0,
  //   tfaceclip_newedges = 0.0,
  //   tfaceclip_final = 0.0,
  //   tfaceclip_deactivate = 0.0,
  //                tcap = 0.0,
  //          tconvertto = 0.0;

  // Convert the polyhedron to faces consisting of edge loops.
  // auto t0 = Timing::currentTime();
  EASY_BLOCK("Convert from Spheral::Polyhedron");
  auto vertices = poly.vertices();             // Note this is a copy!
  vector<Edge> edges;                          // Edges as pairs of vertex indices
  vector<Face> faces;                          // Faces that make up the polyhedron (loops of edges).
  vector<int> vertexMask(vertices.size(), 1);  // Mask to flag active/inactive vertices: -2=>inactive, (-1,0,1)=>(below,in,above) current plane
  vector<int> edgeMask;                        // Mask to flag active/inactive edges: -2=>inactive, -2=>2 points below, 0=>both points in plane, 1=>1 point below/1 above, 2=>both points above
  vector<Vector> faceNormal;                   // Outward pointing unit normals to each face.
  {
    const auto& facets = poly.facets();

    // First pass, build the unique edges in sorted order.
    int iedge, iface = 0;
    for (const auto& facet: facets) {
      const auto& ipoints = facet.ipoints();
      const auto  npoints = ipoints.size();
      for (auto k = 0; k < npoints; ++k) {
        const auto edge = make_edge(ipoints[k], ipoints[(k + 1) % npoints]);
        auto itr = lower_bound(edges.begin(), edges.end(), edge);
        if (itr == edges.end() or *itr != edge) {
          edges.insert(itr, edge);
          edgeMask.push_back(1);
        }
      }
    }

    // Second pass, build the faces as edge loops.
    for (const auto& facet: facets) {
      const auto& ipoints = facet.ipoints();
      const auto  npoints = ipoints.size();
      Face face;
      for (auto k = 0; k < npoints; ++k) {
        const auto edge = make_edge(ipoints[k], ipoints[(k + 1) % npoints]);
        const auto itr = lower_bound(edges.begin(), edges.end(), edge);
        CHECK(itr != edges.end());
        iedge = distance(edges.begin(), itr);
        face.push_back(edge.first == ipoints[k] ? iedge : ~iedge);
      }
      faces.push_back(face);
      faceNormal.push_back(facet.normal());
      ++iface;
    }
  }
  EASY_END_BLOCK;
  // tconvertfrom += Timing::difference(t0, Timing::currentTime());
  
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

    // Check the active vertices against this plane.
    auto above = true;
    auto below = true;
    // t0 = Timing::currentTime();
    EASY_BLOCK("Clip vertices");
    for (auto k = 0; k < vertices.size(); ++k) {
      if (vertexMask[k] != -2) {
        vertexMask[k] = compare(p0, phat, vertices[k]);
        if (vertexMask[k] == 1) {
          below = false;
        } else if (vertexMask[k] == -1) {
          above = false;
        }
      }
    }
    CHECK(not (above and below));
    EASY_END_BLOCK;
    // tvertexmask += Timing::difference(t0, Timing::currentTime());

    // Did we get a simple case?
    if (below) {
      // The current polyhedron is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes either -- we're done.
      edges.clear();
      faces.clear();

    } else if (not above) {

      // This plane passes somewhere through the polyhedron, so we need to clip it.
      // Clip each active edge.
      // t0 = Timing::currentTime();
      EASY_BLOCK("Clip edges");
      vector<int> edgeNodeInPlane(edges.size(), -1);
      vector<int> edgeSign(edges.size(), 1);         // -1=>edge order flipped, 1=>unchanged
      int v0, v1, ivert;
      for (auto iedge = 0; iedge != edges.size(); ++iedge) {
        if (edgeMask[iedge] != -2) {
          tie(v0, v1) = edges[iedge];
          CHECK(vertexMask[v0] != -2 and vertexMask[v1] != -2);

          // Check all the ways this edge could interact with the plane.
          if (vertexMask[v0] == 1 and vertexMask[v1] == 1) {

            // This edge got through unscathed.
            edgeMask[iedge] = 2;

          } else if (vertexMask[v0] == -1 and vertexMask[v1] == 1) {

            // v0 is clipped
            edgeMask[iedge] = 1;
            insertVertex(vertices, ivert, vertexMask, v0, v1, p0, phat);
            edgeNodeInPlane[iedge] = ivert;
            edges[iedge] = make_edge(ivert, v1);
            if (edges[iedge].first == v1) edgeSign[iedge] = -1;

          } else if (vertexMask[v0] == 1 and vertexMask[v1] == -1) {

            // v1 is clipped
            edgeMask[iedge] = 1;
            insertVertex(vertices, ivert, vertexMask, v0, v1, p0, phat);
            edgeNodeInPlane[iedge] = ivert;
            edges[iedge] = make_edge(v0, ivert);
            if (edges[iedge].second == v0) edgeSign[iedge] = -1;

          } else if ((vertexMask[v0] == -1 and vertexMask[v1] <=  0) or
                     (vertexMask[v0] <= 0  and vertexMask[v1] == -1)) {

            // v0 and v1 are clipped
            edgeMask[iedge] = -2;

          } else if (vertexMask[v0] == 0 and vertexMask[v1] == 1) {

            // v0 is in-plane
            edgeMask[iedge] = 1;
            edgeNodeInPlane[iedge] = v0;

          } else if (vertexMask[v0] == 1 and vertexMask[v1] == 0) {

            // v1 is in-plane
            edgeMask[iedge] = 1;
            edgeNodeInPlane[iedge] = v1;

          } else {

            // This edge landed exactly in-plane.
            // In this case we don't add it's nodes to the in-plane node set, since there's already an edge
            // connecting them.
            CHECK2(vertexMask[v0] == 0 and vertexMask[v1] == 0, v0 << " " << v1 << " " << vertexMask[v0] << " " << vertexMask[v1]);
            edgeMask[iedge] = 0;
            // edgeNodesInPlane[iedge].push_back(v0);
            // edgeNodesInPlane[iedge].push_back(v1);

          }
        }
      }
      CHECK(vertexMask.size() == vertices.size());
      CHECK(edgeMask.size() == edges.size());
      EASY_END_BLOCK;
      // tedgemask += Timing::difference(t0, Timing::currentTime());

      // Walk each current face, and reconstruct it according to what happend to it's edges.
      // t0 = Timing::currentTime();
      EASY_BLOCK("Clip faces");
      vector<int> newEdges;                         // Any new edges we create.
      const auto nfaces0 = faces.size();
      for (auto kface = 0; kface < nfaces0; ++kface) {
        // cerr << "Face " << kface << endl;
        // auto t1 = Timing::currentTime();
        EASY_BLOCK("Clip edges in face");
        auto& face = faces[kface];
        Face newface;                  // The new face we're going to build.
        vector<int> faceNodesInPlane;  // Any nodes in this face that are in the clipping plane.

        // Walk the loop of edges in this face.
        for (auto kedge = 0; kedge < face.size(); ++kedge) {
          const auto iedge = posID(face[kedge]);
          if (edgeMask[iedge] != -2) {
            v0 = startVertex(face[kedge], edges);
            v1 = endVertex(face[kedge], edges);
            CHECK(vertexMask[v0] != -2 and vertexMask[v1] != -2);

            if (edgeMask[iedge] == 2) {

              // This edge was unaffected.
              newface.push_back(edgeSign[iedge] == 1 ? face[kedge] : ~face[kedge]);

            } else if (edgeMask[iedge] == 1) {

              // One node of this edge was clipped.
              newface.push_back(edgeSign[iedge] == 1 ? face[kedge] : ~face[kedge]);
              CHECK(edgeNodeInPlane[iedge] >= 0);
              faceNodesInPlane.push_back(edgeNodeInPlane[iedge]);

            } else if (edgeMask[iedge] == 0) {

              // This edge is in the plane.  We don't add it's nodes to the faceNodesInPlane
              // since the edge connecting them is aleady built.
              newface.push_back(edgeSign[iedge] == 1 ? face[kedge] : ~face[kedge]);
              // CHECK(edgeNodesInPlane[iedge].size() == 2);
              // faceNodesInPlane.push_back(edgeNodesInPlane[iedge][0]);
              // faceNodesInPlane.push_back(edgeNodesInPlane[iedge][1]);
            }
          }
        }
        CHECK2(faceNodesInPlane.size() % 2 == 0, faceNodesInPlane.size() << " " << newface.size());   // Gotta come in pairs!
        EASY_END_BLOCK;
        // tfaceclip_edges += Timing::difference(t1, Timing::currentTime());

        // If there's no more than one edge of this face left, the face is gone.
        if (newface.size() <= 1) {
          // t1 = Timing::currentTime();
          EASY_BLOCK("Erase face");
          face.clear();
          EASY_END_BLOCK;
          // tfaceclip_clear += Timing::difference(t1, Timing::currentTime());

        } else {

          // Check for any gaps in the face ring, which will make new edges.
          if (not faceNodesInPlane.empty()) {

            // Get rid of any duplicates.
            // t1 = Timing::currentTime();
            EASY_BLOCK("Erase new node duplicates");
            sort(faceNodesInPlane.begin(), faceNodesInPlane.end());
            vector<int> nodes2kill;
            for (auto k = 0; k < faceNodesInPlane.size() - 1; ++k) {
              if (faceNodesInPlane[k] == faceNodesInPlane[k+1]) {
                nodes2kill.push_back(k);
                nodes2kill.push_back(k + 1);
              }
            }
            removeElements(faceNodesInPlane, nodes2kill);
            EASY_END_BLOCK;
            // tfaceclip_dups += Timing::difference(t1, Timing::currentTime());

            // Sort the hanging nodes in the plane along the line created by the intersection of the
            // t1 = Timing::currentTime();
            EASY_BLOCK("Sort new nodes in line");
            const auto direction = phat.cross(faceNormal[kface]).unitVector();
            sort(faceNodesInPlane.begin(), faceNodesInPlane.end(),
                 [&](const int a, const int b) { return (vertices[a] - p0).dot(direction) < (vertices[b] - p0).dot(direction); });
            EASY_END_BLOCK;
            // tfaceclip_sortline += Timing::difference(t1, Timing::currentTime());

            // Now the ordered pairs of these vertices form new edges to close the rings of the face.
            // t1 = Timing::currentTime();
            EASY_BLOCK("Construct new edges for face");
            for (auto k = 0; k < faceNodesInPlane.size(); k += 2) {
              v0 = faceNodesInPlane[k];
              v1 = faceNodesInPlane[k+1];
              CHECK((vertices[v1] - vertices[v0]).dot(direction) > 0.0);
              const auto edge = make_edge(v0, v1);
              auto itr = find(edges.begin(), edges.end(), edge);
              if (itr == edges.end()) {   // Check if there was an edge embedded in the plane for these nodes already.
                const int iedge = edges.size();
                edges.push_back(edge);
                edgeMask.push_back(1);
                newface.push_back(edges.back().first == v0 ? iedge : ~iedge);
                newEdges.push_back(~newface.back());
              } else {
                const int iedge = distance(edges.begin(), itr);
                CHECK(edgeMask[iedge] == 0);
                newEdges.push_back(edges[iedge].first == v0 ? ~iedge : iedge);
              }
            }
            CHECK(newface.size() >= 3);
            EASY_END_BLOCK;
            // tfaceclip_newedges += Timing::difference(t1, Timing::currentTime());

            // Sort the edges in the face to form contiguous rings.
            // t1 = Timing::currentTime();
            EASY_BLOCK("Construct loops in new face(s)");
            sortInTopologicalRings(newface, edges);

            // Do we need to check for multiple loops in this face?
            if (faceNodesInPlane.size() == 2) {

              // Nope, we have the simple case of just one face.
              face = newface;

            } else {
              // cerr << "Splitting!" << endl;

              // We have clipped a non-convex section of the face, which may have created more than one independent
              // ring of edges -- i.e., we may need to create more than one new face.
              auto kstart = 0;
              while (kstart < newface.size()) {
                v0 = startVertex(newface[kstart], edges);
                auto kend = kstart + 2;  // Gotta be at least three 3 edges per ring.
                while (kend < newface.size() and endVertex(newface[kend], edges) != v0) ++kend;
                CHECK(kend < newface.size());
                Face ring(newface.begin() + kstart, newface.begin() + kend + 1);

                // // BLAGO!
                // {
                //   cerr << "   Ring: ";
                //   for (auto k = 0; k < ring.size(); ++k) {
                //     const auto iedge = posID(ring[k]);
                //     cerr << " ([" << ring[k] << "->" << iedge << "] ";
                //     if (ring[k] < 0) {
                //       cerr << edges[iedge].second << " " << edges[iedge].first << ")";
                //     } else {
                //       cerr << edges[iedge].first << " " << edges[iedge].second << ")";
                //     }
                //   }
                //   cerr << endl;
                // }
                // // BLAGO!
                
                // Verify the new face forms a proper topological ring.
                BEGIN_CONTRACT_SCOPE
                {
                  for (auto kedge = 0; kedge < ring.size(); ++kedge) {
                    const auto kedge1 = (kedge + 1) % ring.size();
                    CHECK(endVertex(ring[kedge], edges) == startVertex(ring[kedge1], edges));
                  }
                }
                END_CONTRACT_SCOPE

                // Which face are we replacing?
                if (kstart == 0) {
                  face = ring;
                } else {
                  faces.push_back(ring);
                  faceNormal.push_back(faceNormal[kface]);
                }

                kstart = kend + 1;
              }
            }
            CHECK(faces.size() == faceNormal.size());
          }
          EASY_END_BLOCK;
          // tfaceclip_final += Timing::difference(t1, Timing::currentTime());
        }

        // Deactivate any clipped vertices.
        // t1 = Timing::currentTime();
        EASY_BLOCK("Mark clipped vertices inactive");
        for (auto k = 0; k < vertexMask.size(); ++k) {
          if (vertexMask[k] == -1) vertexMask[k] = -2;
        }
        EASY_END_BLOCK;
        // tfaceclip_deactivate += Timing::difference(t1, Timing::currentTime());
        // cerr << poly2string(vertices, vertexMask, edges, faces) << endl;
      }
      // tfaceclip += Timing::difference(t0, Timing::currentTime());

      // // BLAGO
      // {
      //   cerr << "newEdges: ";
      //   std::copy(newEdges.begin(), newEdges.end(), std::ostream_iterator<int>(std::cerr, " "));
      //   cerr << endl;
      // }
      // // BLAGO

      // Well, now all the starting faces of the polyhedron have been clipped by the plane, but we need
      // to cap off new edges created by the plane as new faces of the polyhedron.
      // t0 = Timing::currentTime();
      EASY_BLOCK("Cap new face(s) in plane");
      if (not newEdges.empty()) {
        CHECK(newEdges.size() >= 3);

        // We start by ordering the new edges such that they form topological loops.  Note
        // there may be more than one loop if we are cutting a non-convex polyhedron.  We also know the directions
        // we want to walk these loops based on our sign convention for the ID in newEdges.
        sortInTopologicalRings(newEdges, edges);

        // cerr << "Ordered edges for capping: ";
        // for (auto iedge = 0; iedge != newEdges.size(); ++iedge) {
        //   if (newEdges[iedge] >= 0) {
        //     cerr << " (" << edges[newEdges[iedge]].first << " " << edges[newEdges[iedge]].second << ")";
        //   } else {
        //     cerr << " (" << edges[~newEdges[iedge]].second << " " << edges[~newEdges[iedge]].first << ")";
        //   }
        // }
        // cerr << endl;

        auto kstart = 0;
        while (kstart < newEdges.size()) {
          v0 = startVertex(newEdges[kstart], edges);
          auto kend = kstart + 2;  // Gotta be at least three 3 edges per ring.
          while (kend < newEdges.size() and endVertex(newEdges[kend], edges) != v0) ++kend;
          CHECK(kend < newEdges.size());
          Face ring(newEdges.begin() + kstart, newEdges.begin() + kend + 1);
                
          // Verify the new face forms a proper topological ring.
          BEGIN_CONTRACT_SCOPE
          {
            for (auto kedge = 0; kedge < ring.size(); ++kedge) {
              const auto kedge1 = (kedge + 1) % ring.size();
              CHECK(endVertex(ring[kedge], edges) == startVertex(ring[kedge1], edges));
            }
          }
          END_CONTRACT_SCOPE

          // Add the new face.
          faces.push_back(ring);
          faceNormal.push_back(-phat);
          kstart = kend + 1;
        }
      }
      CHECK(faceNormal.size() == faces.size());
      CHECK(vertexMask.size() == vertices.size());
      CHECK(edgeMask.size() == edges.size());
      EASY_END_BLOCK;
      // tcap += Timing::difference(t0, Timing::currentTime());
    }
  }

  // // BLAGO
  // cerr << "----------------------------------------------------------------------" << endl
  //      << "Final polyhedron: " << endl << poly2string(vertices, vertexMask, edges, faces) << endl;
  // // BLAGO

  // If we have an empty polyhedron just finish it here.
  // t0 = Timing::currentTime();
  EASY_BLOCK("Convert to Spheral::Polyhedron");
  if (faces.empty()) {
    poly = GeomPolyhedron();

  } else {

    // Build a list of the vertices that are active, and a map of old->new vertex index.
    CHECK(vertexMask.empty() or 
          (*min_element(vertexMask.begin(), vertexMask.end()) >= -2 and
           *max_element(vertexMask.begin(), vertexMask.end()) <= 1));
    vector<Vector> newvertices;
    newvertices.reserve(vertices.size());
    vector<int> old2new(vertices.size());
    for (auto i = 0; i < vertices.size(); ++i) {
      if (vertexMask[i] >= 0) {
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
    // tconvertto += Timing::difference(t0, Timing::currentTime());
    // cerr << "And the answer is..." << endl
    //      << poly << endl
    //      << "volume = " << poly.volume() << endl;
  }
  EASY_END_BLOCK;

  // // Timing summary.
  // cerr << "Timing summary: read polyhedron: " << tconvertfrom << endl
  //      << "                    vertex mask: " << tvertexmask << endl
  //      << "                      edge mask: " << tedgemask << endl
  //      << "                  face clipping: " << tfaceclip << endl
  //      << "                  --->    edges:   " << tfaceclip_edges << endl
  //      << "                  --->    clear:   " << tfaceclip_clear << endl
  //      << "                  --->     dups:   " << tfaceclip_dups << endl
  //      << "                  ---> sortline:   " << tfaceclip_sortline << endl
  //      << "                  ---> newedges:   " << tfaceclip_newedges << endl
  //      << "                  --->    final:   " << tfaceclip_final << endl
  //      << "                  -> deactivate:   " << tfaceclip_deactivate << endl
  //      << "                   face capping: " << tcap << endl
  //      << "               write polyhedron: " << tconvertto << endl;
}

}
