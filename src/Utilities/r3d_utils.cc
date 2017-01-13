//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using R2D/R3D methods.
//------------------------------------------------------------------------------
#include <algorithm>
#include <set>
#include <iostream>

#include "r3d_utils.hh"

namespace Spheral {

using namespace std;

namespace {   // anonymous namespace

//------------------------------------------------------------------------------
// A class to hold indices making up a planar polygonal face.
// The finalize method shifts the loop to start with the minimum index to make
// each loop unique for comparisons.
//------------------------------------------------------------------------------
struct Face {
  vector<unsigned> indices;
  Face(unsigned i, unsigned j):
    indices({i, j}) {
    CHECK(i != j);
  }
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
Dim<3>::Vector cell_normal(const r3d_poly& celli,
                           const vector<Dim<3>::Vector>& uniqueVerts,
                           const vector<unsigned>& id,
                           unsigned i0,
                           unsigned i1,
                           unsigned i2) {
  return (uniqueVerts[id[i1]] - uniqueVerts[id[i0]]).cross(uniqueVerts[id[i2]] - uniqueVerts[id[i0]]).unitVector();
}

//------------------------------------------------------------------------------
// Walk the R3D vertex connectivity 'til the loop closes.
//------------------------------------------------------------------------------
Face walkR3DFace(const r3d_poly& celli,
                 const vector<Dim<3>::Vector>& uniqueVerts,
                 const vector<unsigned>& id,
                 const unsigned i0,
                 const unsigned i1,
                 const unsigned i2) {

  typedef Dim<3>::Vector Vector;

  Face result(id[i0], id[i1]);

  // Build the normal we check against for the face.
  const Vector fhat = cell_normal(celli, uniqueVerts, id, i0, i1, i2);

  // Walk around the cell topology until we arrive back at the starting vertex.
  unsigned last = i0, next = i1;
  while (id[next] != id[i0]) {
    // Look for the next vertex in the plane.
    CHECK(id[celli.verts[next].pnbrs[0]] == id[last] or
          id[celli.verts[next].pnbrs[1]] == id[last] or
          id[celli.verts[next].pnbrs[2]] == id[last]);
    if (id[celli.verts[next].pnbrs[0]] == id[last]) {                                                                       // Next is either pnbrs[1,2]
      if (abs(abs(cell_normal(celli, uniqueVerts, id, last, next, celli.verts[next].pnbrs[1]).dot(fhat)) - 1.0) < 1.0e-8) { // Check pnbrs[1]
        last = next;
        next = celli.verts[next].pnbrs[1];
      } else {
        last = next;
        next = celli.verts[next].pnbrs[2];
      }
    } else if (id[celli.verts[next].pnbrs[1]] == id[last]) {                                                                // Next is either pnbrs[0,2]
      if (abs(abs(cell_normal(celli, uniqueVerts, id, last, next, celli.verts[next].pnbrs[0]).dot(fhat)) - 1.0) < 1.0e-8) { // Check pnbrs[0]
        last = next;
        next = celli.verts[next].pnbrs[0];
      } else {
        last = next;
        next = celli.verts[next].pnbrs[2];
      }
    } else {                                                                                                                // Next is either pnbrs[0,1]
      if (abs(abs(cell_normal(celli, uniqueVerts, id, last, next, celli.verts[next].pnbrs[0]).dot(fhat)) - 1.0) < 1.0e-8) { // Check pnbrs[0]
        last = next;
        next = celli.verts[next].pnbrs[0];
      } else {
        last = next;
        next = celli.verts[next].pnbrs[1];
      }
    }
    CHECK(id[celli.verts[next].pnbrs[0]] == id[last] or
          id[celli.verts[next].pnbrs[1]] == id[last] or
          id[celli.verts[next].pnbrs[2]] == id[last]);
    CHECK2(id[next] == id[i0] or
           abs(abs(cell_normal(celli, uniqueVerts, id, i0, last, next).dot(fhat)) - 1.0) < 1.0e-8,
           i0 << " " << i1 << " " << i2 << " " << last << " " << next << " " << cell_normal(celli, uniqueVerts, id, i0, last, next) << " " << fhat << " " << cell_normal(celli, uniqueVerts, id, i0, last, next).dot(fhat));
    if (id[next] != id[i0]) result.append(id[next]);
  }

  // Arrange the loop in our standard unique order, and we're done.
  result.finalize();
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
  CHECK(nverts <= R2D_MAX_VERTS);
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
  CHECK(nverts <= R3D_MAX_VERTS);

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

  // if (barf) { // BLAGO
  //   cerr << "Raw verts: " << endl;
  //   for (unsigned j = 0; j != celli.nverts; ++j) {
  //     cerr << " --> " << celli.verts[j].pos.x + ri.x() << " " << celli.verts[j].pos.y + ri.y() << endl;
  //   }
  // } // BLAGO

  // Read out the R2D cell in CCW order.  We have to scan for the positive loop of edges though.
  vector<Vector> verts;
  vector<int> vertcheck(celli.nverts, 0);
  {
    int nextvert, ivert, firstvert;
    double area = -1.0;
    while (area < 0.0) {
      area = 0.0;

      // Find the first unused vertex.
      firstvert = 0;
      while (firstvert != celli.nverts and vertcheck[firstvert] == 1) firstvert++;
      CHECK(firstvert != celli.nverts);

      // Read out the loop of vertices.
      ivert = firstvert;
      nextvert = -1;
      verts.clear();
      while (nextvert != firstvert) {
        verts.push_back(Vector(celli.verts[ivert].pos.x,
                               celli.verts[ivert].pos.y));
        vertcheck[ivert] = 1;
        nextvert = celli.verts[ivert].pnbrs[0];
        // if (barf) cerr << " **> " << (verts.back() + ri) << " " << ivert << " "  << nextvert << endl;
        area += (celli.verts[ivert].pos.x * celli.verts[nextvert].pos.y -
                 celli.verts[ivert].pos.y * celli.verts[nextvert].pos.x);
        ivert = nextvert;
      }
    }
    // if (barf) cerr << " area : " << area << endl;
  }

  // Flag any redundant vertices to not be used.
  vector<int> usevert(verts.size(), 1);
  for (int j = 0; j != verts.size() - 1; ++j) {
    for (int k = j + 1; k != verts.size(); ++k) {
      if (usevert[k] == 1 and (verts[j] - verts[k]).magnitude2() < tol) usevert[k] = 0;
    }
  }

  // Now we can read out the vertices we're actually using and build the return polygon.
  vector<Vector> uniqueVerts;
  vector<vector<unsigned> > facetIndices;
  int k = 0;
  for (int j = 0; j != verts.size(); ++j) {
    if (usevert[j] == 1) {
      uniqueVerts.push_back(verts[j]);
      facetIndices.push_back(vector<unsigned>(2));
      facetIndices.back()[0] = k;
      facetIndices.back()[1] = ++k;
    }
  }
  facetIndices.back()[1] = 0;
  CHECK(uniqueVerts.size() >= 3);

  // // Check the dang things are in CCW order.
  // double area = 0.0;
  // for (int j = 0; j != uniqueVerts.size(); ++j) {
  //   area += ((uniqueVerts[facetIndices[j][0]] - ri).cross(uniqueVerts[facetIndices[j][1]] - ri)).z();
  // }
  // if (area < 0.0) std::reverse(uniqueVerts.begin(), uniqueVerts.end());

  // if (barf) {
  //   cout << " --> " << i << " : ";
  //   std::copy(verts.begin(), verts.end(), std::ostream_iterator<Dim<2>::Vector>(std::cout, " "));
  //   std::cout << endl;
  //   cout << " --> " << i << " : ";
  //   std::copy(uniqueVerts.begin(), uniqueVerts.end(), std::ostream_iterator<Dim<2>::Vector>(std::cout, " "));
  //   std::cout << endl;
  // }

  result = FacetedVolume(uniqueVerts, facetIndices);
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
  typedef std::tuple<unsigned, unsigned, unsigned> face;

  // Is there anything to do?
  result = FacetedVolume();
  if (celli.nverts == 0) return;

  // Build a unique set of vertices.
  vector<Vector> verts;
  vector<unsigned> r3dv2v(celli.nverts);
  unsigned i, j, k;
  for (unsigned i = 0; i != celli.nverts; ++i) {
    const Vector p(celli.verts[i].pos.x,
                   celli.verts[i].pos.y,
                   celli.verts[i].pos.z);
    j = 0;
    while (j < verts.size() and (verts[j] - p).magnitude2() > tol) ++j;
    if (j == verts.size()) verts.push_back(p);
    r3dv2v[i] = j;
    CHECK(j < r3dv2v.size() and r3dv2v[j] < verts.size());
  }

  // Now build the unique faces.
  set<Face> faces;
  for (i = 0; i != celli.nverts; ++i) {
    faces.insert(walkR3DFace(celli, verts, r3dv2v, i, celli.verts[i].pnbrs[0], celli.verts[i].pnbrs[1]));
    faces.insert(walkR3DFace(celli, verts, r3dv2v, i, celli.verts[i].pnbrs[1], celli.verts[i].pnbrs[2]));
    faces.insert(walkR3DFace(celli, verts, r3dv2v, i, celli.verts[i].pnbrs[2], celli.verts[i].pnbrs[0]));
  }
  // sort(faces.begin(), faces.end());
  // const auto lastUniqueFace = unique(faces.begin(), faces.end());

  cerr << "Final faces:" << endl;
  for (auto& face: faces) cerr << face << endl;
  
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
    planes2d[i].d = p.dot(nhat);
  }

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
    planes3d[i].d = p.dot(nhat);
  }

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

}
