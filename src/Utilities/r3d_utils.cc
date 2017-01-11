//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using R2D/R3D methods.
//------------------------------------------------------------------------------
#include "r3d_utils.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Return a Spheral Polygon from an R2D polygon.
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume
r2d_poly_to_polygon(const r2d_poly& celli,
                    const double tol) {

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
  if (celli.nverts == 0) return FacetedVolume();

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

  return FacetedVolume(uniqueVerts, facetIndices);
}

//------------------------------------------------------------------------------
// Return a Spheral Polyhedron from an R3D polyhedron.
//------------------------------------------------------------------------------
Dim<3>::FacetedVolume
r3d_poly_to_polyhedron(const r3d_poly& celli,
                       const double tol) {
  
  using std::vector;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::FacetedVolume::Facet Facet;

  return FacetedVolume();
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

  // Preconditions.
  const auto& vertices = poly.vertices();
  const auto nverts = vertices.size();
  VERIFY2(nverts < R2D_MAX_VERTS,
          "clipFacetedVolume ERROR: input polygon contains " << nverts << " vertices which exceeds max allowed " << R2D_MAX_VERTS);

  // Construct the R2D version of our polygon.
  vector<r2d_rvec2> verts2d(nverts);
  for (unsigned i = 0; i != nverts; ++i) {
    verts2d[i].x = vertices[i].x();
    verts2d[i].y = vertices[i].y();
  }
  r2d_poly poly2d;
  r2d_init_poly(&poly2d, &verts2d[0], nverts);

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
  r2d_real area;
  r2d_reduce(&poly2d, &area, 0);
  const double tol = 1.0e-10 * area;
  return r2d_poly_to_polygon(poly2d, tol);
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

  // Preconditions.
  const auto& vertices = poly.vertices();
  const auto nverts = vertices.size();
  VERIFY2(nverts < R3D_MAX_VERTS,
          "clipFacetedVolume ERROR: input polyhedron contains " << nverts << " vertices which exceeds max allowed " << R3D_MAX_VERTS);

  // Construct the R3D version of our polyhedron.
  const auto& facetVertices = poly.facetVertices();
  const auto nfacets = facetVertices.size();
  vector<r3d_rvec3> verts3d(nverts);
  for (unsigned i = 0; i != nverts; ++i) {
    verts3d[i].x = vertices[i].x();
    verts3d[i].y = vertices[i].y();
    verts3d[i].z = vertices[i].z();
  }
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
  r3d_poly poly3d;
  r3d_init_poly(&poly3d, &verts3d[0], nverts, faceinds, numvertsperface, nfacets);

  // Clean up.
  for (unsigned i = 0; i != nfacets; ++i) {
    delete [] faceinds[i];
  }
  delete [] faceinds;

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
  r3d_real vol;
  r3d_reduce(&poly3d, &vol, 0);
  const double tol = 1.0e-10 * vol;
  return r3d_poly_to_polyhedron(poly3d, tol);
}

}
