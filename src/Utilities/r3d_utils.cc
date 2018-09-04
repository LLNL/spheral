//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using R2D/R3D methods.
//------------------------------------------------------------------------------
#include "r3d_utils.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/sort_permutation.hh"

#include <algorithm>
#include <set>
#include <iostream>
#include <iterator>
using std::vector;
using std::set;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {   // anonymous namespace

//------------------------------------------------------------------------------
// A special comparator to sort r2d planes by distance.
//------------------------------------------------------------------------------
inline
bool compareR2Dplanes(const r2d_plane& lhs, const r2d_plane& rhs) {
  return lhs.d < rhs.d;
}

//------------------------------------------------------------------------------
// A special comparator to sort r3d planes by distance.
//------------------------------------------------------------------------------
inline
bool compareR3Dplanes(const r3d_plane& lhs, const r3d_plane& rhs) {
  return lhs.d < rhs.d;
}

//------------------------------------------------------------------------------
// A class to hold indices making up a planar polygonal face.
// The finalize method shifts the loop to start with the minimum index to make
// each loop unique for comparisons.
//------------------------------------------------------------------------------
struct Face {
  vector<unsigned> indices;
  Face(unsigned i): indices({i}) {}
  void append(const unsigned i) { indices.push_back(i); }
  void finalize() {
    auto minitr = min_element(indices.begin(), indices.end());
    vector<unsigned> tmp(minitr, indices.end());
    tmp.insert(tmp.end(), indices.begin(), minitr);
    indices = tmp;
  }
  bool operator==(const Face& rhs) const { return indices == rhs.indices; }
  bool operator!=(const Face& rhs) const { return indices != rhs.indices; }
  bool operator< (const Face& rhs) const { return indices <  rhs.indices; }
};

std::ostream&
operator<<(std::ostream& os, const Face& face) {
  os << "Face[";
  std::copy(face.indices.begin(), face.indices.end(), std::ostream_iterator<unsigned>(os, " "));
  os << "]";
  return os;
}

//------------------------------------------------------------------------------
// Compute a face normal.
//------------------------------------------------------------------------------
Dim<3>::Vector facet_normal(const vector<Dim<3>::Vector>& verts,
                           unsigned i0,
                           unsigned i1,
                           unsigned i2) {
  return (verts[i1] - verts[i0]).cross(verts[i2] - verts[i0]).unitVector();
}

//------------------------------------------------------------------------------
// Walk the R3D vertex connectivity 'til the loop closes.
//------------------------------------------------------------------------------
Face findFaceRing(const vector<Dim<3>::Vector>& verts,
                  const vector<vector<unsigned> > nghbrs,
                  const GeomPlane<Dim<3> >& facePlane,
                  const unsigned i0,
                  const double tol) {

  typedef Dim<3>::Vector Vector;
  typedef GeomPlane<Dim<3> > Plane;

  // cerr << "----------------------------------------------------------------------" << endl;
  // cerr << i0 << " " << facePlane << endl;

  const Vector& fhat = facePlane.normal();

  // There should be two vertices connected to this one in the plane.  Find the correct
  // one to start the CCW (viewed from the exterior) loop.
  vector<unsigned> otherverts;
  for (auto j: nghbrs[i0]) {
    if (facePlane.minimumDistance(verts[j]) < tol) otherverts.push_back(j);
  }
  CHECK(otherverts.size() == 2);
  unsigned next;
  if (facet_normal(verts, i0, otherverts[0], otherverts[1]).dot(fhat) > 0.0) {
    next = otherverts[0];
  } else {
    next = otherverts[1];
  }

  // cerr << "Stage 1: " << i0 << " " << otherverts[0] << " " << otherverts[1] << " : " << next << endl;

  // Kick off the result.
  Face result(i0);

  // Walk around the cell topology until we arrive back at the starting vertex.
  unsigned last = i0;
  unsigned j, n;
  while (next != i0) {
    result.append(next);

    // Look for the next vertex in the plane.
    j = 0;
    auto itr = nghbrs[next].begin();
    while (itr != nghbrs[next].end() and (*itr == last or facePlane.minimumDistance(verts[*itr]) > tol)) ++itr;
    CHECK(itr != nghbrs[next].end());
    last = next;
    next = *itr;
    // cerr << "  ---> " << last << " " << next << endl;
  }

  // Arrange the loop in our standard unique order, and we're done.
  result.finalize();
  // cerr << "Final ring: " << result << endl;
  // cerr << "          : ";
  // for (const auto i: result.indices) cerr << facePlane.minimumDistance(verts[i]) << " ";
  // cerr << endl;
  // cerr << "          : ";
  // for (const auto i: result.indices) cerr << verts[i] << " ";
  // cerr << endl;
  // cerr << "          : ";
  // cerr << facePlane << " ";
  // cerr << endl;
  return result;
}

}             // anonymous namespace

//------------------------------------------------------------------------------
// Transform a R2D polygon from a Spheral Polygon.
//------------------------------------------------------------------------------
void
polygon_to_r2d_poly(const Dim<2>::FacetedVolume& poly, r2d_poly& result) {

  using std::vector;
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::FacetedVolume::Facet Facet;

  // We use the knowledge here that the polygon vertices are already in CCW order.
  const vector<Vector>& vertices = poly.vertices();
  const unsigned nverts = vertices.size();
  vector<r2d_rvec2> verts2d(nverts);
  for (unsigned i = 0; i != vertices.size(); ++i) {
    verts2d[i].x = vertices[i].x();
    verts2d[i].y = vertices[i].y();
  }
  r2d_init_poly(&result, &verts2d[0], nverts);
}

//------------------------------------------------------------------------------
// Transform a R3D polyhedron from a Spheral Polyhedron.
//------------------------------------------------------------------------------
void
polyhedron_to_r3d_poly(const Dim<3>::FacetedVolume& poly, r3d_poly& result) {

  using std::vector;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::FacetedVolume::Facet Facet;

  const auto& vertices = poly.vertices();
  const auto& facetVertices = poly.facetVertices();
  auto nverts = vertices.size();
  auto nfacets = facetVertices.size();

  // First convert the vertices.
  vector<r3d_rvec3> verts3d(nverts);
  for (unsigned i = 0; i != nverts; ++i) {
    verts3d[i].x = vertices[i].x();
    verts3d[i].y = vertices[i].y();
    verts3d[i].z = vertices[i].z();
  }

  // Now the faces.
  r3d_int numvertsperface[nfacets];
  r3d_int** faceinds = new r3d_int*[nfacets];
  for (unsigned i = 0; i != nfacets; ++i) {
    const unsigned nFacetVerts = facetVertices[i].size();
    numvertsperface[i] = nFacetVerts;
    faceinds[i] = new r3d_int[nFacetVerts];
    for (unsigned j = 0; j != nFacetVerts; ++j) {
      faceinds[i][j] = facetVertices[i][j];
    }
  }
  r3d_init_poly(&result, &verts3d[0], nverts, faceinds, numvertsperface, nfacets);

  // Clean up.
  for (unsigned i = 0; i != nfacets; ++i) {
    delete [] faceinds[i];
  }
  delete [] faceinds;
}

//------------------------------------------------------------------------------
// Return a Spheral Polygon from an R2D polygon.
//------------------------------------------------------------------------------
void
r2d_poly_to_polygon(const r2d_poly& celli,
                    const double tol,
                    Dim<2>::FacetedVolume& result) {

  // Note, R2D leaves lots of degeneracies in the cell points/edges, so we do this in two passes.  First,
  // read all the vertices in CCW order and build a linked list pointing to the next one.  Then we
  // go over these points and remove any degeneracies by updating just the linked list to loop over
  // unique vertices.

  using std::vector;
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::FacetedVolume::Facet Facet;

  // Is there anything to do?
  result = FacetedVolume();
  if (celli.nverts == 0) return;

  // Build the unique set of vertices.
  vector<Vector> verts;
  vector<unsigned> id(celli.nverts);
  unsigned i, j, n;
  for (i = 0; i != celli.nverts; ++i) {
    const Vector p(celli.verts[i].pos.x,
                   celli.verts[i].pos.y);
    j = 0;
    while (j < verts.size() and (verts[j] - p).magnitude2() > tol) ++j;
    if (j == verts.size()) verts.push_back(p);
    id[i] = j;
    CHECK(j < id.size() and id[j] < verts.size());
  }

  // If the input was degenerate, we just kick back an empty polygon.
  const unsigned nunique = verts.size();
  if (nunique < 3) return;

  // Find the centroid.
  // CHECK2(nunique >= 3, "r2d_poly_to_polygon: too few unique nodes: " << nunique);
  const Vector centroid = std::accumulate(verts.begin(), verts.end(), Vector::zero)/nunique;

  // Build the vertex->vertex connectivity.
  vector<vector<unsigned> > nghbrs(nunique);
  for (i = 0; i != celli.nverts; ++i) {
    j = id[i];
    if (id[celli.verts[i].pnbrs[0]] != j) nghbrs[j].push_back(id[celli.verts[i].pnbrs[0]]);
  }
  // for (i = 0; i != celli.nverts; ++i) {
  //   j = id[i];
  //   if (id[celli.verts[i].pnbrs[0]] != j and find(nghbrs[j].begin(), nghbrs[j].end(), id[celli.verts[i].pnbrs[0]]) == nghbrs[j].end()) nghbrs[j].push_back(id[celli.verts[i].pnbrs[0]]);
  // }

  // // BLAGO!
  // {
  //   cerr << "R2D cell:" << endl;
  //   r2d_print(const_cast<r2d_poly*>(&celli));
  //   cerr << "Unique vertices: " << endl;
  //   for (unsigned i = 0; i != nunique; ++i) {
  //     cerr << "                 " << i << "\t" << verts[i] << "\t" << " : ";
  //     std::copy(nghbrs[i].begin(), nghbrs[i].end(), ostream_iterator<unsigned>(cerr, " "));
  //     cerr << endl;
  //   }
  // }
  // // BLAGO!

  // Now read out the final cell in CCW order.  We assume the first loop we find that has postive area
  // (i.e., is in CCW order) is the one we want.
  

  vector<Vector> CCWverts;
  vector<int> vertcheck(nunique, 0);
  {
    unsigned nextvert, ivert, firstvert;
    while (CCWverts.size()< verts.size()) {

      // Find the first unused vertex.
      firstvert = 0;
      while (firstvert != nunique and vertcheck[firstvert] == 1) ++firstvert;
      CHECK(firstvert != nunique);

      // Read out the loop of vertices.
      ivert = firstvert;
      nextvert = nunique;

      while (nextvert != firstvert and CCWverts.size() < verts.size()) {
        CCWverts.push_back(verts[ivert]);
        vertcheck[ivert] = 1;
        nextvert = nghbrs[ivert][0];
        // cerr << " **> " << (CCWverts.back()) << " " << ivert << " "  << nextvert << endl;
        ivert = nextvert;
      }
      // cerr << "Area: " << area << endl;
    }
  }

  const unsigned nCCWverts = CCWverts.size();
  CHECK(nCCWverts >= 3);

  // Now we can build our polygon.
  vector<vector<unsigned> > facetIndices(nCCWverts, vector<unsigned>(2));
  for (unsigned i = 0; i != nCCWverts; ++i) {
    facetIndices[i][0] = i;
    facetIndices[i][1] = (i + 1) % nCCWverts;
  }
  CHECK(facetIndices.back()[1] == 0);
  result = FacetedVolume(CCWverts, facetIndices);
}

//------------------------------------------------------------------------------
// Return a Spheral Polyhedron from an R3D polyhedron.
//------------------------------------------------------------------------------
void
r3d_poly_to_polyhedron(const r3d_poly& celli,
                       const double tol,
                       Dim<3>::FacetedVolume& result) {

  using std::vector;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::FacetedVolume::Facet Facet;
  typedef GeomPlane<Dim<3> > Plane;

  // Is there anything to do?
  result = FacetedVolume();
  if (celli.nverts == 0) return;

  // Build the unique set of vertices.
  vector<Vector> verts;
  vector<unsigned> id(celli.nverts);
  unsigned i, j, k, n;
  for (i = 0; i != celli.nverts; ++i) {
    const Vector p(celli.verts[i].pos.x,
                   celli.verts[i].pos.y,
                   celli.verts[i].pos.z);
    j = 0;
    while (j < verts.size() and (verts[j] - p).magnitude2() > tol) ++j;
    if (j == verts.size()) verts.push_back(p);
    id[i] = j;
    CHECK(j < id.size() and id[j] < verts.size());
  }

  // Build the unique Faces.
  // This bit of code is a version of Joachim's Pouderoux's algorithm.
  set<Face> faces;
  bool skipface;
  unsigned prev, next;
  for (i = 0; i != celli.nverts; ++i) {
    for (k = 0; k != 3; ++k) {
      Face face(id[i]);
      prev = i;
      next = celli.verts[i].pnbrs[k];
      skipface = false;
      do {
        if (id[next] != face.indices.back()) face.append(id[next]);
        j = 0;
        while (j != 3 and celli.verts[next].pnbrs[j] != prev) ++j;
        if (j == 3) {
          skipface = true;
          break;
        }
        prev = next;
        next = celli.verts[prev].pnbrs[(j + 1) % 3];
      } while (prev != i);
      face.indices.pop_back();  // Due to our logic we insert the same index at the end as the start.
      // cerr << "FACE: " << skipface << " : " << face << endl;
      if (not skipface and face.indices.size() >= 3) {
        std::reverse(face.indices.begin(), face.indices.end());  // The above walk seems to get the vertices in CW order, so switch to CCW.
        face.finalize();
        faces.insert(face);
      }
    }
  }

  // // Grab the centroid.
  // const Vector centroid = accumulate(verts.begin(), verts.end(), Vector::zero)/nunique;

  // // Build the vertex->vertex connectivity.
  // const unsigned nunique = verts.size();
  // vector<vector<unsigned> > nghbrs(nunique);
  // vector<vector<Vector> > vertexFaceNormals(nunique);
  // for (i = 0; i != celli.nverts; ++i) {
  //   j = id[i];
  //   if (id[celli.verts[i].pnbrs[0]] != j and find(nghbrs[j].begin(), nghbrs[j].end(), id[celli.verts[i].pnbrs[0]]) == nghbrs[j].end()) nghbrs[j].push_back(id[celli.verts[i].pnbrs[0]]);
  //   if (id[celli.verts[i].pnbrs[1]] != j and find(nghbrs[j].begin(), nghbrs[j].end(), id[celli.verts[i].pnbrs[1]]) == nghbrs[j].end()) nghbrs[j].push_back(id[celli.verts[i].pnbrs[1]]);
  //   if (id[celli.verts[i].pnbrs[2]] != j and find(nghbrs[j].begin(), nghbrs[j].end(), id[celli.verts[i].pnbrs[2]]) == nghbrs[j].end()) nghbrs[j].push_back(id[celli.verts[i].pnbrs[2]]);
  // }

  // // Build the facet normals around each vertex.  This is complicated by the limitations of the r3d conventions for storing topology.
  // for (i = 0; i != nunique; ++i) {
  //   n = nghbrs[i].size();
  //   CHECK(n >= 3);
  //   if (true) { // (n > 3) {
  //     // Oh boy, this point was degenerate so the ordering of neighbors around the vertex is suspect.  We sort them
  //     // in CCW order around point i in a projected plane.  This won't be correct if the local surface is non-convex however.
  //     vector<Vector> npoints;
  //     for (j = 0; j != n; ++j) npoints.push_back(verts[nghbrs[i][j]] - verts[i]);
  //     Plane localplane(npoints);                                                                        // Builds the best-fit local plane
  //     if (localplane.normal().dot(verts[i] - centroid) < 0.0) localplane.normal(-localplane.normal());  // Check normal orientation
  //     CounterClockwiseComparator<Vector, vector<Vector> > CCW(npoints, localplane.point(), localplane.normal());
  //     const auto perm = sort_permutation(npoints, CCW);
  //     apply_permutation_in_place(nghbrs[i], perm);
  //   }
  //   for (j = 0; j != n; ++j) vertexFaceNormals[i].push_back(facet_normal(verts, i, nghbrs[i][j], nghbrs[i][(j+1)%n]));
  // }

  // // BLAGO
  // for (auto& v: verts) cerr << "V---> " << v << endl;
  // for (auto& nbs: nghbrs) {
  //   cerr << "NGB-> ";
  //   for (auto i: nbs) cerr << i << " ";
  //   cerr << endl;
  // }
  // for (auto& normals: vertexFaceNormals) {
  //   cerr << "NRM-> ";
  //   for (auto& n: normals) cerr << n << " ";
  //   cerr << endl;
  // }
  // // BLAGO

  // // Now build the unique faces.
  // set<Face> faces;
  // for (i = 0; i != nunique; ++i) {
  //   for (const auto& nhat: vertexFaceNormals[i]) faces.insert(findFaceRing(verts, nghbrs, Plane(verts[i], nhat), i, tol));
  // }
  // // sort(faces.begin(), faces.end());
  // // const auto lastUniqueFace = unique(faces.begin(), faces.end());

  // cerr << "Final faces:" << endl;
  // for (auto& face: faces) cerr << face << endl;
  
  // Copy the unique faces to a single array.
  vector<vector<unsigned> > faceIndices;
  for (const auto& face: faces) faceIndices.push_back(face.indices);

  // Now we can build the polyhedron.
  result = FacetedVolume(verts, faceIndices);
}

//------------------------------------------------------------------------------
// Clip a polygon by a series of planes
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<2> > >& planes) {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly;

  // Construct the R2D version of our polygon.
  r2d_poly poly2d;
  polygon_to_r2d_poly(poly, poly2d);

  // Now the R2D planes.
  vector<r2d_plane> planes2d(nplanes);
  for (unsigned i = 0; i != nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes2d[i].n.x = nhat.x();
    planes2d[i].n.y = nhat.y();
    planes2d[i].d = -p.dot(nhat);
  }

  // Sort the planes by distance.
  sort(planes2d.begin(), planes2d.end(), compareR2Dplanes);

  // Do the deed.
  r2d_clip(&poly2d, &planes2d[0], nplanes);

  // Copy back to a Spheral polygon.
  FacetedVolume result;
  r2d_real area;
  r2d_reduce(&poly2d, &area, 0);
  const double tol = 1.0e-10 * area;
  r2d_poly_to_polygon(poly2d, tol, result);
  return result;
}

//------------------------------------------------------------------------------
// Clip a polyhedron by a series of planes
//------------------------------------------------------------------------------
Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<3> > >& planes) {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly;

  // Construct the R3D version of our polyhedron.
  r3d_poly poly3d;
  polyhedron_to_r3d_poly(poly, poly3d);

  // Now the R3D planes.
  vector<r3d_plane> planes3d(nplanes);
  for (unsigned i = 0; i != nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes3d[i].n.x = nhat.x();
    planes3d[i].n.y = nhat.y();
    planes3d[i].n.z = nhat.z();
    planes3d[i].d = -p.dot(nhat);
  }

  // Sort the planes by distance.
  sort(planes3d.begin(), planes3d.end(), compareR3Dplanes);

  // Do the deed.
  r3d_clip(&poly3d, &planes3d[0], nplanes);

  // Copy back to a Spheral polyhedron.
  FacetedVolume result;
  r3d_real vol;
  r3d_reduce(&poly3d, &vol, 0);
  const double tol = 1.0e-10 * vol;
  r3d_poly_to_polyhedron(poly3d, tol, result);
  return result;
}

//------------------------------------------------------------------------------
// Volume of the cliped polygon
//------------------------------------------------------------------------------
double clippedVolume(const Dim<2>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<2> > >& planes) {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly.volume();

  // Construct the R2D version of our polygon.
  r2d_poly poly2d;
  polygon_to_r2d_poly(poly, poly2d);

  // Now the R2D planes.
  vector<r2d_plane> planes2d(nplanes);
  for (unsigned i = 0; i != nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes2d[i].n.x = nhat.x();
    planes2d[i].n.y = nhat.y();
    planes2d[i].d = -p.dot(nhat);
  }

  // Sort the planes by distance.
  sort(planes2d.begin(), planes2d.end(), compareR2Dplanes);

  // Do the deed.
  r2d_clip(&poly2d, &planes2d[0], nplanes);

  // Return the volume.
  r2d_real area;
  r2d_reduce(&poly2d, &area, 0);
  return double(area);
}

//------------------------------------------------------------------------------
// Volume of the clipped polyhedron
//------------------------------------------------------------------------------
double clippedVolume(const Dim<3>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<3> > >& planes) {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly.volume();

  // Construct the R3D version of our polyhedron.
  r3d_poly poly3d;
  polyhedron_to_r3d_poly(poly, poly3d);

  // Now the R3D planes.
  vector<r3d_plane> planes3d(nplanes);
  for (unsigned i = 0; i != nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes3d[i].n.x = nhat.x();
    planes3d[i].n.y = nhat.y();
    planes3d[i].n.z = nhat.z();
    planes3d[i].d = -p.dot(nhat);
  }

  // Sort the planes by distance.
  sort(planes3d.begin(), planes3d.end(), compareR3Dplanes);

  // Do the deed.
  r3d_clip(&poly3d, &planes3d[0], nplanes);

  // Return the volume.
  r3d_real vol;
  r3d_reduce(&poly3d, &vol, 0);
  return double(vol);
}

}
