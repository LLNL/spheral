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
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

#include <list>
#include <map>
#include <iostream>
#include <iterator>
#include <algorithm>

// Declare the timers.
extern Timer TIME_PC3d_convertto;
extern Timer TIME_PC3d_convertfrom;
extern Timer TIME_PC3d_copy;
extern Timer TIME_PC3d_moments;
extern Timer TIME_PC3d_clip;
extern Timer TIME_PC3d_planes;
extern Timer TIME_PC3d_checkverts;
extern Timer TIME_PC3d_insertverts;
extern Timer TIME_PC3d_hanging;
extern Timer TIME_PC3d_compress;

namespace PolyClipper {

using namespace std;

namespace {    // anonymous methods

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x > 0 -> sgn0(x) =  1
//    x = 0 -> sgn0(x) =  0
//    x < 0 -> sgn0(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn0(const double x) {
  return (x > 0.0 ?  1.0 :
          x < 0.0 ? -1.0 :
          0.0);
}

//------------------------------------------------------------------------------
// Compare a plane and point (our built-in plane one has some issues).
//------------------------------------------------------------------------------
inline
int compare(const Spheral::Dim<3>::Vector& planePoint,
            const Spheral::Dim<3>::Vector& planeNormal,
            const Spheral::Dim<3>::Vector& point) {
  return sgn0(planeNormal.dot(point - planePoint));
}

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
inline
Spheral::Dim<3>::Vector
segmentPlaneIntersection(const Spheral::Dim<3>::Vector& a,       // line-segment begin
                         const Spheral::Dim<3>::Vector& b,       // line-segment end
                         const Spheral::Dim<3>::Vector& p,       // point in plane
                         const Spheral::Dim<3>::Vector& phat) {  // plane unit normal

  const auto ab = b - a;
  const auto abhat = ab.unitVector();
  CHECK2(std::abs(abhat.dot(phat)) > 0.0, (abhat.dot(phat)) << " " << a << " " << b << " " << abhat << " " << phat);
  const auto s = (p - a).dot(phat)/(abhat.dot(phat));
  CHECK(s >= 0.0 and s <= ab.magnitude());
  const auto result = a + s*abhat;
  // CHECK2(fuzzyEqual((result - p).dot(phat), 0.0, 1.0e-10),
  //        a << " " << b << " " << s << " " << result << " " << (result - p).dot(phat));
  return result;
}

}              // anonymous methods

//------------------------------------------------------------------------------
// Return the vertices ordered in faces.
// Implicitly uses the convention that neighbors for each vertex are arranged
// counter-clockwise viewed from the exterior.
//------------------------------------------------------------------------------
vector<vector<const Vertex3d*>>
extractFaces(const Polyhedron& poly) {

  typedef pair<const Vertex3d*, const Vertex3d*> Edge;
  typedef vector<const Vertex3d*> Face;

  // Prepare the result.
  vector<vector<const Vertex3d*>> result;

  // Numbers of active vertices.
  const auto nactive = count_if(poly.begin(), poly.end(),
                                [](const Vertex3d& x) { return x.comp >= 0; });
  set<Edge> edgesWalked;

  // Walk each vertex in the polyhedron.
  for (const auto& v: poly) {
    if (v.comp >= 0) {

      // Check every (outgoing) edge attached to this vertex.
      for (const auto nptr: v.neighbors) {
        CHECK(nptr->comp >= 0);

        // Has this edge been walked yet?
        if (edgesWalked.find(make_pair(&v, nptr)) == edgesWalked.end()) {

          // Follow around the face represented by this edge until we get back
          // to our starting vertex.
          Face face;
          face.push_back(&v);
          auto vnext = nptr;
          auto vprev = &v;
          while (vnext != &v) {
            face.push_back(vnext);
            CHECK(edgesWalked.find(make_pair(vprev, vnext)) == edgesWalked.end());
            edgesWalked.insert(make_pair(vprev, vnext));
            auto itr = find(vnext->neighbors.begin(), vnext->neighbors.end(), vprev);
            CHECK(itr != vnext->neighbors.end());
            vprev = vnext;
            if (itr == vnext->neighbors.begin()) {
              vnext = vnext->neighbors.back();
            } else {
              vnext = *(itr - 1);
            }
          }
          CHECK(face.size() >= 3);
          result.push_back(face);

        }
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    // Every pair should have been walked twice, once in each direction.
    for (const auto& v: poly) {
      for (const auto nptr: v.neighbors) {
        CHECK(edgesWalked.find(make_pair(&v, nptr)) != edgesWalked.end());
        CHECK(edgesWalked.find(make_pair(nptr, &v)) != edgesWalked.end());
      }
    }
  }
  END_CONTRACT_SCOPE

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polyhedron.
//------------------------------------------------------------------------------
std::string
polyhedron2string(const Polyhedron& poly) {
  ostringstream s;
  s << "[";

  // Get the vertices in face ordering.
  const vector<vector<const Vertex3d*>> faces;

  // Now output the face vertex coordinates.
  for (const auto& face: faces) {
    s << "[";
    for (const auto vptr: face) {
      s << " " << vptr->position;
    }
    s << "]\n ";
  }
  s << "]";

  return s.str();
}

//------------------------------------------------------------------------------
// Convert Spheral::GeomPolyhedron -> PolyClipper::Polyhedron.
//------------------------------------------------------------------------------
void convertToPolyhedron(Polyhedron& polyhedron,
                      const Spheral::Dim<3>::FacetedVolume& Spheral_polyhedron) {
  TIME_PC3d_convertto.start();

  const auto& vertPositions = Spheral_polyhedron.vertices();
  const auto& facets = Spheral_polyhedron.facets();
  const auto  nverts = vertPositions.size();
  const auto  nfacets = facets.size();

  // Build the PolyClipper Vertex3d's, but without connectivity yet.
  Polyhedron result;
  vector<Vertex3d*> id2vert;
  for (auto k = 0; k < nverts; ++k) {
    result.push_back(Vertex3d(vertPositions[k], 1));
    id2vert.push_back(&result.back());
  }
  CHECK(id2vert.size() == nverts);

  // Note all the edges associated with each vertex.
  vector<vector<pair<int, int>>> vertexPairs(nverts);
  int v0, vprev, vnext;
  for (const auto& facet: facets) {
    const auto& ipoints = facet.ipoints();
    const int   n = ipoints.size();
    for (int k = 0; k < n; ++k) {
      v0 = ipoints[k];
      vprev = ipoints[(k - 1) % n];
      vnext = ipoints[(k + 1) % n];
      vertexPairs[v0].push_back(make_pair(vnext, vprev));   // Gotta reverse!
    }
  }

  // Now we can build vertex->vertex connectivity.
  for (auto k = 0; k < nverts; ++k) {

    // Sort the edges associated with this vertex.
    const auto n = vertexPairs[k].size();
    CHECK(n >= 3);
    for (auto i = 0; i < n - 1; ++i) {
      auto j = i + 1;
      while (j < n and vertexPairs[k][i].second != vertexPairs[k][j].first) ++j;
      CHECK(j < n);
      swap(vertexPairs[k][i + 1], vertexPairs[k][j]);
    }

    // Now that they're in order, create the neighbors.
    for (auto i = 0; i < n; ++i) {
      id2vert[k]->neighbors.push_back(id2vert[vertexPairs[k][i].first]);
    }
    CHECK(id2vert[k]->neighbors.size() == vertexPairs[k].size());
  }

  CHECK(polyhedron.size() == nverts);
  TIME_PC3d_convertto.stop();
}

// //------------------------------------------------------------------------------
// // Convert PolyClipper::Polyhedron -> Spheral::GeomPolyhedron.
// //------------------------------------------------------------------------------
// void convertFromPolyhedron(Spheral::Dim<3>::FacetedVolume& Spheral_polyhedron,
//                         const Polyhedron& polyhedron) {
//   TIME_PC3d_convertfrom.start();

//   // Useful types.
//   typedef Spheral::Dim<3>::FacetedVolume FacetedVolume;
//   typedef Spheral::Dim<3>::Vector Vector;

//   if (polyhedron.empty()) {

//     Spheral_polyhedron = FacetedVolume();

//   } else {

//     // Numbers of vertices.
//     const auto nactive = count_if(polyhedron.begin(), polyhedron.end(),
//                                   [](const Vertex3d& x) { return x.comp >= 0; });
//     set<const Vertex3d*> usedVertices;

//     // Go until we hit all the active vertices.
//     vector<Vector> coords(nactive);
//     vector<vector<unsigned>> facets(nactive, vector<unsigned>(2));
//     auto k = 0, loopStart = 0;
//     while (usedVertices.size() < nactive) {

//       // Look for the first active unused vertex.
//       auto vstart = polyhedron.begin();
//       while (vstart != polyhedron.end() and
//              (vstart->comp < 0 or usedVertices.find(&(*vstart)) != usedVertices.end())) vstart++;
//       CHECK(vstart != polyhedron.end());
//       auto vnext = &(*vstart);

//       // Read out this loop.
//       auto force = true;
//       while (force or vnext != &(*vstart)) {
//         CHECK(k < nactive);
//         coords[k] = vnext->position;
//         facets[k][0] = k;
//         facets[k][1] = k + 1;
//         ++k;
//         force = false;
//         usedVertices.insert(vnext);
//         vnext = vnext->neighbors.second;
//       }
//       facets[k-1][1] = loopStart;
//       loopStart = k;
//     }
//     CHECK(k == nactive);

//     Spheral_polyhedron = FacetedVolume(coords, facets);

//   }
//   TIME_PC3d_convertfrom.stop();
// }

// //------------------------------------------------------------------------------
// // Copy a PolyClipper::Polyhedron.
// //------------------------------------------------------------------------------
// void copyPolyhedron(Polyhedron& polyhedron,
//                  const Polyhedron& polyhedron0) {
//   TIME_PC3d_copy.start();
//   polyhedron.clear();
//   if (not polyhedron0.empty()) {
//     std::map<const Vertex3d*, Vertex3d*> ptrMap;
//     for (auto& v: polyhedron0) {
//       polyhedron.push_back(v);
//       ptrMap[&v] = &polyhedron.back();
//     }
//     for (auto& v: polyhedron) {
//       v.neighbors.first  = ptrMap[v.neighbors.first];
//       v.neighbors.second = ptrMap[v.neighbors.second];
//     }
//   }
//   TIME_PC3d_copy.stop();
// }

// //------------------------------------------------------------------------------
// // Compute the zeroth and first moment of a Polyhedron.
// //------------------------------------------------------------------------------
// void moments(double& zerothMoment, Spheral::Dim<3>::Vector& firstMoment,
//              const Polyhedron& polyhedron) {
//   TIME_PC3d_moments.start();

//   // Useful types.
//   typedef Spheral::Dim<3>::Vector Vector;

//   // Clear the result for accumulation.
//   zerothMoment = 0.0;
//   firstMoment = Vector::zero;

//   // Walk the polyhedron, and add up our results triangle by triangle.
//   if (not polyhedron.empty()) {
//     const auto nverts = polyhedron.size();
//     auto vfirst = &polyhedron.front();
//     auto vptr = vfirst->neighbors.second;
//     Vertex3d* vnext;
//     for (auto k = 0; k < nverts; ++k) {
//       vnext = vptr->neighbors.second;
//       const auto triA = (vptr->position - vfirst->position).cross(vnext->position - vfirst->position).z();
//       zerothMoment += triA;
//       firstMoment += triA * (vptr->position + vnext->position);
//       vptr = vnext;
//     }
//     CHECK(zerothMoment != 0.0);
//     firstMoment = firstMoment/(3.0*zerothMoment) + vfirst->position;
//     zerothMoment *= 0.5;
//   }
//   TIME_PC3d_moments.stop();
// }

// //------------------------------------------------------------------------------
// // Clip a polyhedron by planes.
// //------------------------------------------------------------------------------
// void clipPolyhedron(Polyhedron& polyhedron,
//                  const std::vector<Spheral::GeomPlane<Spheral::Dim<3>>>& planes) {
//   TIME_PC3d_clip.start();

//   // Useful types.
//   typedef Spheral::Dim<3>::Vector Vector;

//   // Loop over the planes.
//   TIME_PC3d_planes.start();
//   auto kplane = 0;
//   const auto nplanes = planes.size();
//   while (kplane < nplanes and not polyhedron.empty()) {
//     const auto& plane = planes[kplane++];
//     const auto& p0 = plane.point();
//     const auto& phat = plane.normal();
//     // cerr << "Clip plane: " << p0 << " " << phat << endl;

//     // Check the current set of vertices against this plane.
//     TIME_PC3d_checkverts.start();
//     auto above = true;
//     auto below = true;
//     for (auto& v: polyhedron) {
//       v.comp = compare(p0, phat, v.position);
//       if (v.comp >= 0) {
//         below = false;
//       } else if (v.comp == -1) {
//         above = false;
//       }
//     }
//     CHECK(not (above and below));
//     TIME_PC3d_checkverts.stop();

//     // Did we get a simple case?
//     if (below) {
//       // The polyhedron is entirely below the clip plane, and is therefore entirely removed.
//       // No need to check any more clipping planes -- we're done.
//       polyhedron.clear();

//     } else if (not above) {

//       // This plane passes through the polyhedron.
//       // Insert any new vertices.
//       TIME_PC3d_insertverts.start();
//       vector<Vertex3d*> hangingVertices;
//       Vertex3d *vprev, *vnext;
//       for (auto& v: polyhedron) {
//         std::tie(vprev, vnext) = v.neighbors;

//         if ((v.comp)*(vnext->comp) == -1) {
//           // This pair straddles the plane and creates a new vertex.
//           polyhedron.push_back(Vertex3d(segmentPlaneIntersection(v.position,
//                                                               vnext->position,
//                                                               p0,
//                                                               phat),
//                                      2));         // 2 indicates new vertex
//           polyhedron.back().neighbors.first = &v;
//           polyhedron.back().neighbors.second = vnext;
//           v.neighbors.second = &polyhedron.back();
//           vnext->neighbors.first = &polyhedron.back();
//           hangingVertices.push_back(&polyhedron.back());
//           // cerr << " --> Inserting new vertex @ " << polyhedron.back().position << endl;

//         } else if (v.comp == 0 and 
//                    (vprev->comp == -1 xor vnext->comp == -1)) {
//           // This vertex is exactly in-plane, but has exactly one neighbor edge that will be entirely clipped.
//           // No new vertex, but vptr will be hanging.
//           hangingVertices.push_back(&v);
//           // cerr << " --> Hanging vertex @ " << v.position << endl;

//         }
//       }
//       TIME_PC3d_insertverts.stop();

//       // For each hanging vertex, link to the neighbors that survive the clipping.
//       // If there are more than two hanging vertices, we've clipped a non-convex face and need to check
//       // how to hook up each section, possibly resulting in new faces.
//       TIME_PC3d_hanging.start();
//       CHECK(hangingVertices.size() % 2 == 0);
//       if (true) { //(hangingVertices.size() > 2) {

//         // Yep, more than one new edge here.
//         const Vector direction(phat.y(), -phat.x());
//         sort(hangingVertices.begin(), hangingVertices.end(), 
//              [&](const Vertex3d* a, const Vertex3d* b) { return (a->position - p0).dot(direction) < (b->position - p0).dot(direction); });

//         // Now the ordered pairs of these new vertices form the new edges.
//         int v0, v1, kedge0, kedge1, iedge;
//         Vertex3d *v1ptr, *v2ptr;
//         for (auto k = 0; k < hangingVertices.size(); k += 2) {
//           v1ptr = hangingVertices[k];
//           v2ptr = hangingVertices[k + 1];
//           v1ptr->neighbors.second = v2ptr;
//           v2ptr->neighbors.first = v1ptr;
//         }

//       } else {

//         // Just hook across the vertices and we're done.
//         for (auto vptr: hangingVertices) {
//           std::tie(vprev, vnext) = vptr->neighbors;
//           CHECK(vptr->comp == 0 or vptr->comp == 2);
//           CHECK(vprev->comp == -1 xor vnext->comp == -1);

//           if (vprev->comp == -1) {
//             // We have to search backwards.
//             while (vprev->comp == -1) {
//               vprev = vprev->neighbors.first;
//             }
//             CHECK(vprev != vptr);
//             vptr->neighbors.first = vprev;

//           } else {
//             // We have to search forward.
//             while (vnext->comp == -1) {
//               vnext = vnext->neighbors.second;
//             }
//             CHECK(vnext != vptr);
//             vptr->neighbors.second = vnext;

//           }
//         }

//       }

//       // Remove the clipped vertices, compressing the polyhedron.
//       TIME_PC3d_compress.start();
//       for (auto vitr = polyhedron.begin(); vitr != polyhedron.end();) {
//         if (vitr->comp < 0) {
//           vitr = polyhedron.erase(vitr);
//         } else {
//           ++vitr;
//         }
//       }
//       TIME_PC3d_compress.stop();

//       // cerr << "After compression: " << polyhedron2string(polyhedron) << endl;

//       // Is the polyhedron hedrone?
//       if (polyhedron.size() < 3) polyhedron.clear();
//       TIME_PC3d_hanging.stop();
//     }
//   }
//   TIME_PC3d_planes.stop();
// }

}
